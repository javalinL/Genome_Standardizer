# -*- coding: utf-8 -*-
"""
backend.py — Async pipeline wrapper for gstd TUI.

Provides a non-blocking interface over gstd's synchronous plugin functions.
All blocking calls are executed in a ThreadPoolExecutor; log output is
intercepted via a thread-safe Queue and streamed back to the TUI in real time.

──────────────────────────────────────────────────────────────────────────────
Message Protocol
──────────────────────────────────────────────────────────────────────────────
GstdBackend.run() is an async generator that yields dicts:

  Progress update (emitted before each pipeline step):
    { 'type': 'progress', 'step': int, 'total': int, 'label': str }

  Log line (intercepted from the gstd logger):
    { 'type': 'log', 'level': str, 'message': str }

  Successful completion:
    { 'type': 'result',
      'genes':   int,   # number of genes standardized
      'cds':     int,   # number of CDS sequences extracted
      'pep':     int,   # number of PEP sequences extracted
      'elapsed': float  # wall-clock seconds
    }
    { 'type': 'done' }

  Unrecoverable error (pipeline halted):
    { 'type': 'error', 'message': str }

──────────────────────────────────────────────────────────────────────────────
Usage (inside a Textual Worker or async task)
──────────────────────────────────────────────────────────────────────────────
    from genome_standardizer.tui.backend import GstdBackend, PipelineArgs

    backend = GstdBackend()
    args = PipelineArgs(
        gff_list   = ["annotation.gff3.gz"],
        fasta_list = ["genome.fasta.gz"],
        prefix     = "Oryz_sati",
    )
    async for msg in backend.run(args):
        if msg['type'] == 'progress':
            progress_bar.update(msg['step'], msg['total'])
            step_label.update(msg['label'])
        elif msg['type'] == 'log':
            rich_log.write(msg['message'])
        elif msg['type'] == 'result':
            status_bar.update(f"Done: {msg['genes']} genes, {msg['cds']} CDS")
        elif msg['type'] == 'error':
            rich_log.write(f"[ERROR] {msg['message']}")
        elif msg['type'] == 'done':
            break
"""

import asyncio
import logging
import os
import queue as stdlib_queue
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from typing import AsyncGenerator, List, Optional

from genome_standardizer.plugins import (
    annot_processor,
    cleanup_plugin,
    exporter_plugin,
    parser_plugin,
    seq_extractor,
)
from genome_standardizer.main import build_robust_alias_map

# Total pipeline steps (used by TUI to size the progress bar)
TOTAL_STEPS = 6


# ── Pipeline arguments dataclass ─────────────────────────────────────────────

@dataclass
class PipelineArgs:
    """
    Clean, argparse-independent representation of all gstd parameters.

    The TUI builds one of these from its widget values and passes it to
    GstdBackend.run().  Mirrors the fields of the argparse.Namespace
    produced by main.parse_args() so that plugin calls are identical.
    """
    # Required
    gff_list:    List[str]
    fasta_list:  List[str]
    prefix:      str

    # Optional — match CLI defaults
    prefix_list: List[str] = field(default_factory=list)
    step:        int        = 10
    longest:     bool       = False
    keep:        bool       = False
    save_log:    bool       = False
    no_repair:   bool       = False
    pep_qc:      bool       = False
    bed6:        bool       = False
    keep_source: bool       = False
    low_mem:     bool       = False


# ── Thread-safe logging handler ───────────────────────────────────────────────

class _ThreadQueueHandler(logging.Handler):
    """
    Routes log records from a background thread into a stdlib.queue.Queue
    so that the async event loop can drain them without blocking.
    """

    def __init__(self, q: stdlib_queue.Queue) -> None:
        super().__init__()
        self._q = q

    def emit(self, record: logging.LogRecord) -> None:
        try:
            self._q.put_nowait({
                "level":   record.levelname,
                "message": self.format(record).rstrip(),
            })
        except Exception:
            self.handleError(record)


# ── Backend ───────────────────────────────────────────────────────────────────

class GstdBackend:
    """
    Async wrapper around the gstd plugin pipeline.

    One instance can be shared across the lifetime of the TUI application.
    Each call to run() creates its own thread-pool future and log queue,
    so concurrent runs are isolated (though the TUI should only allow one
    at a time in practice).
    """

    def __init__(self, max_workers: int = 1) -> None:
        self._executor = ThreadPoolExecutor(
            max_workers=max_workers,
            thread_name_prefix="gstd_worker",
        )

    # ── Public async generator ────────────────────────────────────────────────

    async def run(self, args: PipelineArgs) -> AsyncGenerator[dict, None]:
        """
        Execute the full gstd pipeline asynchronously.

        Yields message dicts (see module docstring for the full protocol).
        The caller must consume all messages until 'done' or 'error' is
        received.
        """
        loop       = asyncio.get_running_loop()
        log_q      = stdlib_queue.Queue()
        start_time = time.time()

        # ── Wire up logging interception ──────────────────────────────────────
        gstd_logger = logging.getLogger("gstd")
        interceptor = _ThreadQueueHandler(log_q)
        interceptor.setFormatter(logging.Formatter("%(message)s"))
        gstd_logger.addHandler(interceptor)

        file_handler: Optional[logging.FileHandler] = None
        if args.save_log:
            work_dir     = os.path.dirname(os.path.abspath(args.gff_list[0]))
            log_path     = os.path.join(work_dir, f"{args.prefix}_run.log")
            file_handler = logging.FileHandler(log_path, encoding="utf-8")
            file_handler.setFormatter(
                logging.Formatter(
                    "%(asctime)s [%(levelname)s] %(message)s",
                    datefmt="%Y-%m-%d %H:%M:%S",
                )
            )
            gstd_logger.addHandler(file_handler)

        # ── Shared state across steps ─────────────────────────────────────────
        # Each step writes its output into these; the next step reads from them.
        # Using a single-element list as a mutable cell (avoids nonlocal).
        _genome_seqs   = [None]
        _raw_genes     = [None]
        _global_alias  = [None]   # dict
        _processed     = [None]
        _id_mapping    = [None]
        _cds_records   = [None]
        _pep_records   = [None]
        _n_genes       = [0]
        _n_cds         = [0]
        _work_dir      = [None]

        try:
            # ═══════════════════════════════════════════════════════════════
            # STEP 1 — Parse inputs
            # ═══════════════════════════════════════════════════════════════
            yield {
                "type":  "progress",
                "step":  1,
                "total": TOTAL_STEPS,
                "label": "[1/6] Loading and parsing inputs…",
            }

            def _step1() -> None:
                seqs, genes = parser_plugin.parse_inputs(
                    args.gff_list,
                    args.fasta_list,
                    args.prefix_list,
                    gstd_logger,
                    no_repair=args.no_repair,
                    low_mem=args.low_mem,
                )
                _genome_seqs[0] = seqs
                _raw_genes[0]   = genes
                _work_dir[0]    = os.path.dirname(
                    os.path.abspath(args.gff_list[0])
                )

            error = await self._run_step(loop, log_q, _step1)
            async for msg in self._drain(log_q):
                yield msg
            if error:
                yield {"type": "error", "message": str(error)}
                return

            yield {
                "type":    "log",
                "level":   "INFO",
                "message": (
                    f"         Loaded {len(list(_genome_seqs[0].keys()))} scaffolds "
                    f"and parsed {len(_raw_genes[0])} genes."
                ),
            }

            # ═══════════════════════════════════════════════════════════════
            # STEP 2 — Smart Alias Sniffer
            # ═══════════════════════════════════════════════════════════════
            yield {
                "type":  "progress",
                "step":  2,
                "total": TOTAL_STEPS,
                "label": "[2/6] Running Smart Alias Sniffer…",
            }

            def _step2() -> None:
                alias_map   = {}
                added       = 0
                genome_seqs = _genome_seqs[0]

                for fasta_path in args.fasta_list:
                    amap, conflicts = build_robust_alias_map(fasta_path)
                    alias_map.update(amap)
                    if conflicts:
                        gstd_logger.warning(
                            f"         [Warning] {len(conflicts)} conflicting aliases "
                            f"in {os.path.basename(fasta_path)}."
                        )

                for alias, real_id in alias_map.items():
                    if real_id in genome_seqs and alias not in genome_seqs:
                        if args.low_mem:
                            genome_seqs.add_alias(alias, real_id)
                        else:
                            genome_seqs[alias] = genome_seqs[real_id]
                        added += 1

                _global_alias[0] = alias_map
                if added > 0:
                    gstd_logger.info(
                        f"         [Auto-Heal] Established {added} smart mapping links."
                    )
                else:
                    gstd_logger.info(
                        "         [Status OK] No sequence name inconsistencies detected."
                    )

            error = await self._run_step(loop, log_q, _step2)
            async for msg in self._drain(log_q):
                yield msg
            if error:
                yield {"type": "error", "message": str(error)}
                return

            # ═══════════════════════════════════════════════════════════════
            # STEP 3 — Standardize & rename annotations
            # ═══════════════════════════════════════════════════════════════
            yield {
                "type":  "progress",
                "step":  3,
                "total": TOTAL_STEPS,
                "label": "[3/6] Structuring & renaming annotations…",
            }

            def _step3() -> None:
                processed, mapping = annot_processor.standardize_and_rename(
                    _raw_genes[0],
                    prefix=args.prefix,
                    step=args.step,
                    longest_only=args.longest,
                    keep_source=args.keep_source,
                )
                _processed[0]  = processed
                _id_mapping[0] = mapping
                _n_genes[0]    = len(processed)

            error = await self._run_step(loop, log_q, _step3)
            async for msg in self._drain(log_q):
                yield msg
            if error:
                yield {"type": "error", "message": str(error)}
                return

            yield {
                "type":    "log",
                "level":   "INFO",
                "message": (
                    f"         Successfully standardized {_n_genes[0]} genes."
                ),
            }

            # ═══════════════════════════════════════════════════════════════
            # STEP 4 — Extract CDS / PEP sequences
            # ═══════════════════════════════════════════════════════════════
            yield {
                "type":  "progress",
                "step":  4,
                "total": TOTAL_STEPS,
                "label": "[4/6] Extracting CDS / PEP sequences…",
            }

            def _step4() -> None:
                cds, pep = seq_extractor.extract_all(
                    _genome_seqs[0],
                    _processed[0],
                    gstd_logger,
                    pep_qc=args.pep_qc,
                )
                _cds_records[0] = cds
                _pep_records[0] = pep
                _n_cds[0]       = len(cds)

            error = await self._run_step(loop, log_q, _step4)
            async for msg in self._drain(log_q):
                yield msg
            if error:
                yield {"type": "error", "message": str(error)}
                return

            yield {
                "type":    "log",
                "level":   "INFO",
                "message": (
                    f"         Extracted {_n_cds[0]} CDS and "
                    f"{len(_pep_records[0])} PEP sequences."
                ),
            }

            # ═══════════════════════════════════════════════════════════════
            # STEP 5 — Export standardized files
            # ═══════════════════════════════════════════════════════════════
            yield {
                "type":  "progress",
                "step":  5,
                "total": TOTAL_STEPS,
                "label": "[5/6] Exporting standardized files…",
            }

            def _step5() -> None:
                exporter_plugin.export_all(
                    _processed[0],
                    _cds_records[0],
                    _pep_records[0],
                    _id_mapping[0],
                    prefix=args.prefix,
                    work_dir=_work_dir[0],
                    logger=gstd_logger,
                    bed6=args.bed6,
                    keep_source=args.keep_source,
                )

            error = await self._run_step(loop, log_q, _step5)
            async for msg in self._drain(log_q):
                yield msg
            if error:
                yield {"type": "error", "message": str(error)}
                return

            # ═══════════════════════════════════════════════════════════════
            # STEP 6 — Compress originals (optional)
            # ═══════════════════════════════════════════════════════════════
            yield {
                "type":  "progress",
                "step":  6,
                "total": TOTAL_STEPS,
                "label": "[6/6] Compressing original inputs…" if not args.keep
                         else "[6/6] Skipping compression (--keep).",
            }

            if not args.keep:
                def _step6() -> None:
                    cleanup_plugin.compress_inputs(
                        args.gff_list + args.fasta_list, gstd_logger
                    )
                error = await self._run_step(loop, log_q, _step6)
                async for msg in self._drain(log_q):
                    yield msg
                if error:
                    yield {"type": "error", "message": str(error)}
                    return
            else:
                yield {
                    "type":    "log",
                    "level":   "INFO",
                    "message": "         [Step 6] Skipped compression (--keep applied).",
                }

            # ═══════════════════════════════════════════════════════════════
            # SUCCESS
            # ═══════════════════════════════════════════════════════════════
            elapsed = time.time() - start_time
            yield {
                "type":    "result",
                "genes":   _n_genes[0],
                "cds":     _n_cds[0],
                "pep":     len(_pep_records[0]) if _pep_records[0] else 0,
                "elapsed": elapsed,
            }
            yield {"type": "done"}

        finally:
            # Always clean up the log interceptor, even on error
            gstd_logger.removeHandler(interceptor)
            if file_handler:
                gstd_logger.removeHandler(file_handler)
                file_handler.close()

    # ── Private helpers ───────────────────────────────────────────────────────

    async def _run_step(
        self,
        loop: asyncio.AbstractEventLoop,
        log_q: stdlib_queue.Queue,
        fn,
        *fn_args,
    ) -> Optional[Exception]:
        """
        Run *fn* in the thread-pool executor.

        Polls every 50 ms while the future is pending, draining *log_q*
        after each sleep (handled by the caller via _drain).

        Returns:
            None on success, or the Exception on failure.
        """
        # Wrap fn so exceptions are captured rather than propagated from
        # the executor future (which would make it hard to yield 'error').
        exc_cell = [None]

        def _safe_run() -> None:
            try:
                fn(*fn_args)
            except Exception as exc:            # noqa: BLE001
                exc_cell[0] = exc

        future = loop.run_in_executor(self._executor, _safe_run)

        while not future.done():
            await asyncio.sleep(0.05)

        await future   # propagate any unexpected executor-level error
        return exc_cell[0]

    @staticmethod
    async def _drain(log_q: stdlib_queue.Queue) -> AsyncGenerator[dict, None]:
        """Drain all pending log messages from *log_q* as 'log' message dicts."""
        while not log_q.empty():
            try:
                record = log_q.get_nowait()
                yield {"type": "log", **record}
            except stdlib_queue.Empty:
                break
