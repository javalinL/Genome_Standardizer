# !/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import argparse
import time
import os
import re
import gzip
import logging

try:
    from genome_standardizer.plugins import parser_plugin
    from genome_standardizer.plugins import annot_processor
    from genome_standardizer.plugins import seq_extractor
    from genome_standardizer.plugins import exporter_plugin
    from genome_standardizer.plugins import cleanup_plugin
except ImportError as e:
    sys.exit(f"[FATAL] Failed to load modules: {e}")


def setup_logger(prefix, work_dir, save_log):
    logger = logging.getLogger("gstd")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    console_fmt = logging.Formatter('%(message)s')
    ch = logging.StreamHandler()
    ch.setFormatter(console_fmt)
    logger.addHandler(ch)

    if save_log:
        file_fmt = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
        log_file = os.path.join(work_dir, f"{prefix}_run.log")
        fh = logging.FileHandler(log_file, encoding='utf-8')
        fh.setFormatter(file_fmt)
        logger.addHandler(fh)

    return logger


def build_robust_alias_map(fasta_file):
    alias_map = {}
    conflict_set = set()
    opener = gzip.open if str(fasta_file).endswith('.gz') else open

    with opener(fasta_file, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                header_parts = line[1:].strip().split(maxsplit=1)
                real_id = header_parts[0]
                alias_map[real_id] = real_id

                if len(header_parts) > 1:
                    description = header_parts[1]
                    generated_aliases = set()

                    chr_match = re.search(r'chromosome\s+([A-Za-z0-9_]+)', description, re.IGNORECASE)
                    if chr_match:
                        val = chr_match.group(1)
                        generated_aliases.update([
                            f"{val}", f"Chr{val}", f"chr{val}",
                            f"Chr{val}.2", f"chr{val}.2",
                            f"Chr{val}_hap1", f"Chr{val}_hap2"
                        ])

                    scaf_match = re.search(r'scaffold[_\s]+([A-Za-z0-9_.]+)', description, re.IGNORECASE)
                    if scaf_match:
                        val = scaf_match.group(1)
                        generated_aliases.update([f"scaffold_{val}", f"Scaffold_{val}"])

                    for alias in generated_aliases:
                        if alias in conflict_set:
                            continue
                        if alias in alias_map:
                            if alias_map[alias] != real_id:
                                del alias_map[alias]
                                conflict_set.add(alias)
                        else:
                            alias_map[alias] = real_id
    return alias_map, conflict_set


def parse_args():
    desc_text = """
Program:  gstd (Genome Standardizer)
Version:  2.4.0 (Stable Release)
Summary:  A robust pipeline for standardizing plant genomes, resolving annotation 
          inconsistencies, and extracting CDS/PEP sequences cleanly.
"""
    parser = argparse.ArgumentParser(
        description=desc_text,
        formatter_class=argparse.RawTextHelpFormatter,
        usage="gstd <GFF> <FASTA> <PREFIX> [OPTIONS]"
    )

    parser.add_argument("gff", help="Input GFF3/GTF file(s). Use comma to separate multiple subgenomes.")
    parser.add_argument("fasta", help="Input FASTA file(s). Use comma to separate multiple subgenomes.")
    parser.add_argument("prefix", help="Target output prefix (e.g., Oryz_sati).")

    group = parser.add_argument_group("Optional Arguments")
    group.add_argument("--step", type=int, default=10, help="Step size for gene numbering (default: 10)")
    group.add_argument("--add-prefix", type=str, default="", help="Isolate subgenomes safely (e.g., SubA_,SubB_)")
    group.add_argument("--longest", action="store_true", help="Keep only the longest transcript per gene")
    group.add_argument("--keep", action="store_true", help="Keep original input files (skip compression)")
    group.add_argument("--save-log", action="store_true", help="Generate a detailed .log file in the output directory")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    args.gff_list = args.gff.split(',')
    args.fasta_list = args.fasta.split(',')

    if len(args.gff_list) != len(args.fasta_list):
        sys.exit(f"[FATAL] GFF count ({len(args.gff_list)}) does not match FASTA count ({len(args.fasta_list)}).")

    args.prefix_list = args.add_prefix.split(',') if args.add_prefix else []
    if args.prefix_list and len(args.prefix_list) != len(args.gff_list):
        sys.exit(f"[FATAL] Isolated prefix count does not match input files count.")

    return args


def main():
    args = parse_args()
    start_time = time.time()
    work_dir = os.path.dirname(os.path.abspath(args.gff_list[0]))

    logger = setup_logger(args.prefix, work_dir, args.save_log)

    logger.info("=" * 60)
    logger.info(" Genome Standardization Pipeline Initialized")
    logger.info(f" Target Prefix : {args.prefix}")
    logger.info(f" Working Dir   : {work_dir}")
    if args.save_log:
        logger.info(f" Log File      : {os.path.join(work_dir, args.prefix + '_run.log')}")
    if len(args.gff_list) > 1:
        mode_text = f"Polyploid Fusion ({len(args.gff_list)} parts)"
        if args.prefix_list:
            mode_text += f" [Isolated with: {', '.join(args.prefix_list)}]"
        logger.info(f" Mode          : {mode_text}")
    if args.longest:
        logger.info(" Transcript    : Longest Only")
    else:
        logger.info(" Transcript    : All Isoforms")
    logger.info("=" * 60)

    try:
        logger.info("\n[Step 1] Loading and Parsing Inputs...")
        genome_seqs, raw_genes = parser_plugin.parse_inputs(args.gff_list, args.fasta_list, args.prefix_list, logger)
        logger.info(f"         Loaded {len(genome_seqs)} scaffolds and parsed {len(raw_genes)} genes.")

        logger.info("\n[Step 1.5] Running Smart Alias Sniffer...")
        global_alias_map = {}
        for fasta in args.fasta_list:
            amap, conflicts = build_robust_alias_map(fasta)
            global_alias_map.update(amap)
            if conflicts:
                logger.warning(
                    f"         [Warning] Intercepted {len(conflicts)} conflicting aliases in {os.path.basename(fasta)}.")

        added_aliases = 0
        for alias, real_id in global_alias_map.items():
            if real_id in genome_seqs and alias not in genome_seqs:
                genome_seqs[alias] = genome_seqs[real_id]
                added_aliases += 1

        if added_aliases > 0:
            logger.info(f"         [Auto-Heal] Established {added_aliases} smart mapping links for memory sequences.")
        else:
            logger.info(f"         [Status OK] No sequence name inconsistencies detected.")

        logger.info("\n[Step 2] Processing Annotations (Structuring & Renaming)...")
        processed_genes, id_mapping = annot_processor.standardize_and_rename(
            raw_genes, prefix=args.prefix, step=args.step, longest_only=args.longest
        )
        logger.info(f"         Successfully standardized {len(processed_genes)} genes.")

        logger.info("\n[Step 3] Extracting Sequences (CDS, PEP)...")
        cds_records, pep_records = seq_extractor.extract_all(genome_seqs, processed_genes, logger)
        logger.info(f"         Extracted {len(cds_records)} CDS and {len(pep_records)} PEP sequences.")

        logger.info("\n[Step 4] Exporting Standardized Files...")
        exporter_plugin.export_all(
            processed_genes, cds_records, pep_records, id_mapping,
            prefix=args.prefix, work_dir=work_dir, logger=logger
        )

        if not args.keep:
            logger.info("\n[Step 5] Compressing Original Inputs...")
            cleanup_plugin.compress_inputs(args.gff_list + args.fasta_list, logger)
        else:
            logger.info("\n[Step 5] Skipped compression (--keep applied).")

        elapsed = time.time() - start_time
        logger.info(f"\n[Done] Standardization completed successfully in {elapsed:.2f} seconds.")

    except Exception as e:
        logger.error(f"\n[FATAL ERROR] Pipeline halted unexpectedly:\n  {e}")
        logger.error("Suggestion: Check input format or use --save-log for detailed diagnostics.")
        sys.exit(1)


if __name__ == "__main__":
    main()