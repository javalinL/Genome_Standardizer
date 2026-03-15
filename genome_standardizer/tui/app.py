# -*- coding: utf-8 -*-
"""
app.py — gstd Terminal User Interface (Fully Implemented by Gemini).
"""

from __future__ import annotations

import asyncio
from typing import Optional
import os

from textual.app import App, ComposeResult
from textual.binding import Binding
from textual.containers import Grid, Horizontal, Vertical, VerticalScroll
from textual.widgets import (
    Button, Collapsible, DirectoryTree, Footer, Header, Input,
    Label, ProgressBar, RichLog, Select, Switch, TabbedContent, TabPane, Markdown
)
from textual.screen import ModalScreen

from genome_standardizer.tui.backend import GstdBackend, PipelineArgs
from genome_standardizer.tui.job_generator import HpcJobParams, generate_script
from genome_standardizer.tui.i18n import I18n, SUPPORTED_LANGS

DEFAULT_LANG = "en"
APP_VERSION = "4.0.0"


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  MODAL SCREENS (File Browser & Help)                                     ║
# ╚═══════════════════════════════════════════════════════════════════════════╝

class FileBrowserModal(ModalScreen[str]):
    def compose(self) -> ComposeResult:
        with Vertical(id="file_browser_container"):
            yield Label(self.app.i18n.t("file_browser.title"), classes="section-title", id="fb_title")
            yield Label(self.app.i18n.t("file_browser.hint"), id="fb_hint")
            self._dir_tree = DirectoryTree(os.path.expanduser("~"), id="file_tree")
            yield self._dir_tree
            with Horizontal(classes="modal-buttons"):
                yield Button(self.app.i18n.t("buttons.cancel"), id="cancel_btn", variant="error")
                yield Button(self.app.i18n.t("buttons.confirm"), id="confirm_btn", variant="success")

    def on_button_pressed(self, event: Button.Pressed) -> None:
        if event.button.id == "cancel_btn":
            self.dismiss(None)
        elif event.button.id == "confirm_btn":
            if self._dir_tree.cursor_node and self._dir_tree.cursor_node.data:
                self.dismiss(self._dir_tree.cursor_node.data.path)
            else:
                self.dismiss(None)

    def on_directory_tree_file_selected(self, event: DirectoryTree.FileSelected) -> None:
        # User pressed Space/Enter on a file
        self.dismiss(event.path)


class HelpModal(ModalScreen):
    def compose(self) -> ComposeResult:
        with Vertical(id="help_container"):
            yield Label(self.app.i18n.t("help.title"), classes="section-title", id="help_title")
            with VerticalScroll():
                yield Markdown(self.app.i18n.t("help.body"), id="help_body")
            with Horizontal(classes="modal-buttons"):
                yield Button(self.app.i18n.t("buttons.close"), id="close_btn", variant="primary")

    def on_button_pressed(self, event: Button.Pressed) -> None:
        if event.button.id == "close_btn":
            self.dismiss()


# ╔═══════════════════════════════════════════════════════════════════════════╗
# ║  MAIN APP (Fully Assembled)                                              ║
# ╚═══════════════════════════════════════════════════════════════════════════╝

class GstdTuiApp(App):
    CSS_PATH = "app.tcss"
    TITLE = f"gstd — Genome Standardizer  v{APP_VERSION}"

    BINDINGS = [
        Binding("ctrl+q", "quit", "Quit", show=True, priority=True),
        Binding("ctrl+l", "toggle_lang", "Language", show=True),
        Binding("ctrl+r", "run", "Run", show=True),
        Binding("ctrl+s", "gen_script", "Gen Script", show=True),
        Binding("f1", "help", "Help", show=True),
    ]

    def __init__(self) -> None:
        super().__init__()
        self.i18n: I18n = I18n(DEFAULT_LANG)
        self.backend: GstdBackend = GstdBackend()
        self._running: bool = False
        self._last_script: str = ""

    def compose(self) -> ComposeResult:
        yield Header()

        with TabbedContent(id="main_tabs"):
            # ── 1. RUN TAB ──────────────────────────────────────────────────
            with TabPane(self.i18n.t("tabs.run"), id="tab_run"):
                yield Label(self.i18n.t("input.section_title"), classes="section-title", id="input_sec_title")
                with Horizontal(classes="field-row"):
                    yield Label(self.i18n.t("input.gff_label"), id="gff_label", classes="field-label")
                    yield Input(placeholder=self.i18n.t("input.gff_placeholder"), id="gff_input")
                    yield Button(self.i18n.t("input.browse_btn"), id="gff_browse_btn")

                with Horizontal(classes="field-row"):
                    yield Label(self.i18n.t("input.fasta_label"), id="fasta_label", classes="field-label")
                    yield Input(placeholder=self.i18n.t("input.fasta_placeholder"), id="fasta_input")
                    yield Button(self.i18n.t("input.browse_btn"), id="fasta_browse_btn")

                with Horizontal(classes="field-row"):
                    yield Label(self.i18n.t("input.prefix_label"), id="prefix_label", classes="field-label")
                    yield Input(placeholder=self.i18n.t("input.prefix_placeholder"), id="prefix_input")

                yield Label(self.i18n.t("input.multi_file_hint"), id="multi_file_hint", classes="text-muted")

                # Basic Options
                yield Label(self.i18n.t("basic_options.section_title"), classes="section-title", id="basic_sec_title")
                with Horizontal(classes="field-row"):
                    yield Label(self.i18n.t("basic_options.step_label"), id="step_label", classes="field-label")
                    yield Input(value="10", id="step_input", type="integer", classes="step-input")
                with Horizontal(classes="toggles-container"):
                    yield Switch(value=False, id="longest_switch")
                    yield Label(self.i18n.t("basic_options.longest_label"), id="longest_label", classes="switch-label")
                    yield Switch(value=False, id="keep_switch")
                    yield Label(self.i18n.t("basic_options.keep_label"), id="keep_label", classes="switch-label")
                    yield Switch(value=False, id="save_log_switch")
                    yield Label(self.i18n.t("basic_options.save_log_label"), id="save_log_label",
                                classes="switch-label")

                # Advanced Collapsible
                with Collapsible(title=self.i18n.t("advanced_options.section_title"), id="advanced_collapsible"):
                    with Horizontal(classes="field-row"):
                        yield Label(self.i18n.t("advanced_options.add_prefix_label"), id="add_prefix_label",
                                    classes="field-label")
                        yield Input(placeholder=self.i18n.t("advanced_options.add_prefix_placeholder"),
                                    id="add_prefix_input")
                    with Horizontal(classes="toggles-container"):
                        yield Switch(value=False, id="no_repair_switch")
                        yield Label(self.i18n.t("advanced_options.no_repair_label"), id="no_repair_label",
                                    classes="switch-label")
                        yield Switch(value=False, id="pep_qc_switch")
                        yield Label(self.i18n.t("advanced_options.pep_qc_label"), id="pep_qc_label",
                                    classes="switch-label")
                        yield Switch(value=False, id="bed6_switch")
                        yield Label(self.i18n.t("advanced_options.bed6_label"), id="bed6_label", classes="switch-label")
                    with Horizontal(classes="toggles-container"):
                        yield Switch(value=False, id="keep_source_switch")
                        yield Label(self.i18n.t("advanced_options.keep_source_label"), id="keep_source_label",
                                    classes="switch-label")
                        yield Switch(value=False, id="low_mem_switch")
                        yield Label(self.i18n.t("advanced_options.low_mem_label"), id="low_mem_label",
                                    classes="switch-label")

                with Horizontal(classes="action-buttons"):
                    yield Button(self.i18n.t("buttons.run_local"), id="run_btn")
                    yield Button(self.i18n.t("buttons.gen_script"), id="gen_script_btn")

            # ── 2. LOG TAB ──────────────────────────────────────────────────
            with TabPane(self.i18n.t("tabs.log"), id="tab_log"):
                yield Label(self.i18n.t("log.idle"), id="progress_label")
                yield ProgressBar(total=6, show_eta=False, id="progress_bar")
                yield RichLog(id="rich_log", markup=True, highlight=True)
                with Horizontal(classes="action-buttons"):
                    yield Button(self.i18n.t("buttons.clear_log"), id="clear_log_btn")
                    yield Button(self.i18n.t("buttons.copy_log"), id="copy_log_btn")

            # ── 3. JOB SCRIPT TAB ───────────────────────────────────────────
            with TabPane(self.i18n.t("tabs.job_script"), id="tab_job"):
                yield Label(self.i18n.t("job_script.intro"), id="job_intro", classes="text-muted")

                with Horizontal(classes="input-row"):
                    with Vertical():
                        yield Label(self.i18n.t("job_script.scheduler_label"), id="scheduler_label")
                        yield Select([("SLURM", "slurm"), ("PBS", "pbs")], value="slurm", id="scheduler_select")
                    with Vertical():
                        yield Label(self.i18n.t("job_script.job_name_label"), id="job_name_label")
                        yield Input(value="gstd_run", id="job_name_input")
                    with Vertical():
                        yield Label(self.i18n.t("job_script.partition_label"), id="partition_label")
                        yield Input(value="normal", id="partition_input")

                with Horizontal(classes="input-row"):
                    with Vertical():
                        yield Label(self.i18n.t("job_script.nodes_label"), id="nodes_label")
                        yield Input(value="1", id="nodes_input", type="integer")
                    with Vertical():
                        yield Label(self.i18n.t("job_script.cpus_label"), id="cpus_label")
                        yield Input(value="8", id="cpus_input", type="integer")
                    with Vertical():
                        yield Label(self.i18n.t("job_script.mem_label"), id="mem_label")
                        yield Input(value="32", id="mem_input", type="integer")
                    with Vertical():
                        yield Label(self.i18n.t("job_script.walltime_label"), id="walltime_label")
                        yield Input(value="12:00:00", id="walltime_input")

                with Horizontal(classes="input-row"):
                    with Vertical():
                        yield Label(self.i18n.t("job_script.conda_env_label"), id="conda_env_label")
                        yield Input(value="core_dev", id="conda_env_input")
                    with Vertical():
                        yield Label(self.i18n.t("job_script.work_dir_label"), id="work_dir_label")
                        yield Input(value="/scratch/user/project", id="work_dir_input")
                    with Vertical():
                        yield Label(self.i18n.t("job_script.output_log_label"), id="output_log_label")
                        yield Input(value="gstd_%j.log", id="output_log_input")

                with Horizontal(classes="action-buttons"):
                    yield Button(self.i18n.t("job_script.gen_btn"), id="gen_job_btn")
                    yield Button(self.i18n.t("buttons.copy_script"), id="copy_script_btn", disabled=True)
                    yield Button(self.i18n.t("buttons.save_script"), id="save_script_btn", disabled=True)

                yield Label(self.i18n.t("job_script.preview_title"), classes="section-title", id="preview_title")
                yield RichLog(id="script_preview", highlight=True)

        yield Footer()

    # ── MOUNT & TOOLTIPS ────────────────────────────────────────────────────
    def on_mount(self) -> None:
        self._apply_tooltips()

    def _apply_tooltips(self) -> None:
        """Assign tooltips post-construction (tooltip= is not a valid ctor kwarg)."""
        self.query_one("#step_label", Label).tooltip       = self.i18n.t("basic_options.step_tooltip")
        self.query_one("#longest_switch", Switch).tooltip  = self.i18n.t("basic_options.longest_tooltip")
        self.query_one("#keep_switch", Switch).tooltip     = self.i18n.t("basic_options.keep_tooltip")
        self.query_one("#save_log_switch", Switch).tooltip = self.i18n.t("basic_options.save_log_tooltip")
        self.query_one("#add_prefix_label", Label).tooltip   = self.i18n.t("advanced_options.add_prefix_tooltip")
        self.query_one("#no_repair_switch", Switch).tooltip  = self.i18n.t("advanced_options.no_repair_tooltip")
        self.query_one("#pep_qc_switch", Switch).tooltip     = self.i18n.t("advanced_options.pep_qc_tooltip")
        self.query_one("#bed6_switch", Switch).tooltip       = self.i18n.t("advanced_options.bed6_tooltip")
        self.query_one("#keep_source_switch", Switch).tooltip = self.i18n.t("advanced_options.keep_source_tooltip")
        self.query_one("#low_mem_switch", Switch).tooltip    = self.i18n.t("advanced_options.low_mem_tooltip")

    # ── EVENTS & LOGIC ───────────────────────────────────────────────────────
    def on_button_pressed(self, event: Button.Pressed) -> None:
        if event.button.id == "gff_browse_btn":
            def check_gff(path: str | None) -> None:
                if path: self.query_one("#gff_input", Input).value = path

            self.push_screen(FileBrowserModal(), check_gff)

        elif event.button.id == "fasta_browse_btn":
            def check_fa(path: str | None) -> None:
                if path: self.query_one("#fasta_input", Input).value = path

            self.push_screen(FileBrowserModal(), check_fa)

        elif event.button.id == "run_btn":
            self.action_run()

        elif event.button.id == "gen_script_btn":
            self.query_one(TabbedContent).active = "tab_job"

        elif event.button.id == "gen_job_btn":
            self.action_gen_script()

        elif event.button.id == "clear_log_btn":
            self.query_one("#rich_log", RichLog).clear()
            self.query_one("#progress_bar", ProgressBar).progress = 0
            self.query_one("#progress_label", Label).update(self.i18n.t("log.idle"))

        elif event.button.id == "copy_log_btn":
            self._copy_to_clipboard(self._richlog_text("#rich_log"), "log_copied")

        elif event.button.id == "copy_script_btn":
            self._copy_to_clipboard(self._last_script, "script_copied")

        elif event.button.id == "save_script_btn":
            self._save_script()

    def _richlog_text(self, widget_id: str) -> str:
        """Extract plain text from a RichLog widget."""
        log = self.query_one(widget_id, RichLog)
        return "\n".join("".join(seg.text for seg in strip._segments) for strip in log.lines)

    def _copy_to_clipboard(self, text: str, notify_key: str) -> None:
        try:
            import pyperclip
            pyperclip.copy(text)
            self.notify(self.i18n.t(f"notify.{notify_key}"))
        except Exception as e:
            self.notify(self.i18n.t("notify.clipboard_unavailable", error=e), severity="warning")

    def _save_script(self) -> None:
        job_name = self.query_one("#job_name_input", Input).value or "gstd_run"
        work_dir = self.query_one("#work_dir_input", Input).value or "."
        path = os.path.join(work_dir, f"{job_name}.sh")
        try:
            with open(path, "w") as f:
                f.write(self._last_script)
            self.notify(self.i18n.t("notify.script_saved", path=path))
        except Exception as e:
            self.notify(self.i18n.t("notify.save_failed", error=e), severity="error")

    def action_toggle_lang(self) -> None:
        current_idx = SUPPORTED_LANGS.index(self.i18n.lang)
        next_lang = SUPPORTED_LANGS[(current_idx + 1) % len(SUPPORTED_LANGS)]
        self.i18n.set_lang(next_lang)
        self._refresh_labels()
        self.notify(self.i18n.t("notify.lang_switched", lang=next_lang))

    def action_run(self) -> None:
        err = self._validate_inputs()
        if err:
            self.notify(err, severity="error")
            return
        args = self._collect_pipeline_args()
        self.query_one(TabbedContent).active = "tab_log"
        self._start_pipeline(args)

    def action_gen_script(self) -> None:
        err = self._validate_inputs()
        if err:
            self.notify(err, severity="error")
            return
        job = self._collect_job_params()
        self._last_script = generate_script(job)
        preview = self.query_one("#script_preview", RichLog)
        preview.clear()
        preview.write(self._last_script)
        self.query_one("#copy_script_btn", Button).disabled = False
        self.query_one("#save_script_btn", Button).disabled = False
        self.notify(self.i18n.t("job_script.generated_hint"))

    def action_help(self) -> None:
        self.push_screen(HelpModal())

    def _validate_inputs(self) -> Optional[str]:
        gff = self.query_one("#gff_input", Input).value.strip()
        fa = self.query_one("#fasta_input", Input).value.strip()
        prefix = self.query_one("#prefix_input", Input).value.strip()
        if not gff: return self.i18n.t("validation.gff_required")
        if not fa: return self.i18n.t("validation.fasta_required")
        if not prefix: return self.i18n.t("validation.prefix_required")
        return None

    def _collect_pipeline_args(self) -> PipelineArgs:
        return PipelineArgs(
            gff_list=[x.strip() for x in self.query_one("#gff_input", Input).value.split(",") if x.strip()],
            fasta_list=[x.strip() for x in self.query_one("#fasta_input", Input).value.split(",") if x.strip()],
            prefix=self.query_one("#prefix_input", Input).value.strip(),
            prefix_list=[x.strip() for x in self.query_one("#add_prefix_input", Input).value.split(",") if x.strip()],
            step=int(self.query_one("#step_input", Input).value or 10),
            longest=self.query_one("#longest_switch", Switch).value,
            keep=self.query_one("#keep_switch", Switch).value,
            save_log=self.query_one("#save_log_switch", Switch).value,
            no_repair=self.query_one("#no_repair_switch", Switch).value,
            pep_qc=self.query_one("#pep_qc_switch", Switch).value,
            bed6=self.query_one("#bed6_switch", Switch).value,
            keep_source=self.query_one("#keep_source_switch", Switch).value,
            low_mem=self.query_one("#low_mem_switch", Switch).value,
        )

    def _collect_job_params(self) -> HpcJobParams:
        args = self._collect_pipeline_args()
        return HpcJobParams(
            scheduler=self.query_one("#scheduler_select", Select).value,
            job_name=self.query_one("#job_name_input", Input).value,
            partition=self.query_one("#partition_input", Input).value,
            nodes=int(self.query_one("#nodes_input", Input).value or 1),
            cpus=int(self.query_one("#cpus_input", Input).value or 8),
            mem_gb=int(self.query_one("#mem_input", Input).value or 32),
            walltime=self.query_one("#walltime_input", Input).value,
            conda_env=self.query_one("#conda_env_input", Input).value,
            work_dir=self.query_one("#work_dir_input", Input).value,
            output_log=self.query_one("#output_log_input", Input).value,
            gff=",".join(args.gff_list),
            fasta=",".join(args.fasta_list),
            prefix=args.prefix,
            extra_flags=" ".join([
                "--longest" if args.longest else "",
                "--keep" if args.keep else "",
                "--save-log" if args.save_log else "",
                "--no-repair" if args.no_repair else "",
                "--pep-qc" if args.pep_qc else "",
                "--bed6" if args.bed6 else "",
                "--low-mem" if args.low_mem else "",
                f"--add-prefix {','.join(args.prefix_list)}" if args.prefix_list else ""
            ]).strip()
        )

    def _refresh_labels(self) -> None:
        # ── Tab titles ────────────────────────────────────────────────────────
        tabs = self.query_one("#main_tabs", TabbedContent)
        tabs.get_tab("tab_run").label = self.i18n.t("tabs.run")
        tabs.get_tab("tab_log").label = self.i18n.t("tabs.log")
        tabs.get_tab("tab_job").label = self.i18n.t("tabs.job_script")
        self._apply_tooltips()

        # ── Run tab — Input section ───────────────────────────────────────────
        self.query_one("#input_sec_title", Label).update(self.i18n.t("input.section_title"))
        self.query_one("#gff_label",   Label).update(self.i18n.t("input.gff_label"))
        self.query_one("#fasta_label", Label).update(self.i18n.t("input.fasta_label"))
        self.query_one("#prefix_label",Label).update(self.i18n.t("input.prefix_label"))
        self.query_one("#multi_file_hint", Label).update(self.i18n.t("input.multi_file_hint"))
        for btn in self.query("Button#gff_browse_btn, Button#fasta_browse_btn"):
            btn.label = self.i18n.t("input.browse_btn")
        self.query_one("#gff_input",   Input).placeholder = self.i18n.t("input.gff_placeholder")
        self.query_one("#fasta_input", Input).placeholder = self.i18n.t("input.fasta_placeholder")
        self.query_one("#prefix_input",Input).placeholder = self.i18n.t("input.prefix_placeholder")

        # ── Run tab — Basic Options ───────────────────────────────────────────
        self.query_one("#basic_sec_title",  Label).update(self.i18n.t("basic_options.section_title"))
        self.query_one("#step_label",       Label).update(self.i18n.t("basic_options.step_label"))
        self.query_one("#longest_label",    Label).update(self.i18n.t("basic_options.longest_label"))
        self.query_one("#keep_label",       Label).update(self.i18n.t("basic_options.keep_label"))
        self.query_one("#save_log_label",   Label).update(self.i18n.t("basic_options.save_log_label"))

        # ── Run tab — Advanced Options ────────────────────────────────────────
        self.query_one("#advanced_collapsible", Collapsible).title = self.i18n.t("advanced_options.section_title")
        self.query_one("#add_prefix_label",  Label).update(self.i18n.t("advanced_options.add_prefix_label"))
        self.query_one("#no_repair_label",   Label).update(self.i18n.t("advanced_options.no_repair_label"))
        self.query_one("#pep_qc_label",      Label).update(self.i18n.t("advanced_options.pep_qc_label"))
        self.query_one("#bed6_label",        Label).update(self.i18n.t("advanced_options.bed6_label"))
        self.query_one("#keep_source_label", Label).update(self.i18n.t("advanced_options.keep_source_label"))
        self.query_one("#low_mem_label",     Label).update(self.i18n.t("advanced_options.low_mem_label"))
        self.query_one("#add_prefix_input",  Input).placeholder = self.i18n.t("advanced_options.add_prefix_placeholder")

        # ── Run tab — Buttons ────────────────────────────────────────────────
        self.query_one("#run_btn",       Button).label = self.i18n.t("buttons.run_local")
        self.query_one("#gen_script_btn",Button).label = self.i18n.t("buttons.gen_script")

        # ── Log tab ───────────────────────────────────────────────────────────
        if not self._running:
            self.query_one("#progress_label", Label).update(self.i18n.t("log.idle"))
        self.query_one("#clear_log_btn", Button).label = self.i18n.t("buttons.clear_log")
        self.query_one("#copy_log_btn",  Button).label = self.i18n.t("buttons.copy_log")

        # ── Job Script tab ────────────────────────────────────────────────────
        self.query_one("#job_intro",        Label).update(self.i18n.t("job_script.intro"))
        self.query_one("#scheduler_label",  Label).update(self.i18n.t("job_script.scheduler_label"))
        self.query_one("#job_name_label",   Label).update(self.i18n.t("job_script.job_name_label"))
        self.query_one("#partition_label",  Label).update(self.i18n.t("job_script.partition_label"))
        self.query_one("#nodes_label",      Label).update(self.i18n.t("job_script.nodes_label"))
        self.query_one("#cpus_label",       Label).update(self.i18n.t("job_script.cpus_label"))
        self.query_one("#mem_label",        Label).update(self.i18n.t("job_script.mem_label"))
        self.query_one("#walltime_label",   Label).update(self.i18n.t("job_script.walltime_label"))
        self.query_one("#conda_env_label",  Label).update(self.i18n.t("job_script.conda_env_label"))
        self.query_one("#work_dir_label",   Label).update(self.i18n.t("job_script.work_dir_label"))
        self.query_one("#output_log_label", Label).update(self.i18n.t("job_script.output_log_label"))
        self.query_one("#preview_title",    Label).update(self.i18n.t("job_script.preview_title"))
        self.query_one("#gen_job_btn",      Button).label = self.i18n.t("job_script.gen_btn")
        self.query_one("#copy_script_btn",  Button).label = self.i18n.t("buttons.copy_script")
        self.query_one("#save_script_btn",  Button).label = self.i18n.t("buttons.save_script")

    def _start_pipeline(self, args: PipelineArgs) -> None:
        if self._running:
            self.notify(self.i18n.t("log.pipeline_running"), severity="warning")
            return
        self._running = True
        self.run_worker(self._execute_run(args), exclusive=True)

    async def _execute_run(self, args: PipelineArgs) -> None:
        try:
            async for msg in self.backend.run(args):
                msg_type = msg.get("type")
                if msg_type == "progress":
                    try:
                        self.query_one("#progress_bar", ProgressBar).update(progress=msg["step"], total=msg["total"])
                        self.query_one("#progress_label", Label).update(msg["label"])
                    except Exception:
                        pass
                elif msg_type == "log":
                    try:
                        clean = msg["message"].rstrip()
                        if clean:
                            self.query_one("#rich_log", RichLog).write(clean)
                    except Exception:
                        pass
                elif msg_type == "result":
                    self._on_run_complete(msg)
                elif msg_type == "error":
                    self._on_run_error(msg["message"])
                    return
                elif msg_type == "done":
                    break
        finally:
            self._running = False

    def _on_run_complete(self, result: dict) -> None:
        summary = self.i18n.t("log.result_summary", genes=result.get("genes", 0), cds=result.get("cds", 0),
                              pep=result.get("pep", 0))
        done_msg = self.i18n.t("log.done", elapsed=result.get("elapsed", 0.0))
        self.notify(done_msg, severity="information")
        try:
            log = self.query_one("#rich_log", RichLog)
            log.write(done_msg)
            log.write(summary)
        except Exception:
            pass

    def _on_run_error(self, message: str) -> None:
        err_msg = self.i18n.t("log.error", message=message)
        self.notify(err_msg, severity="error")
        try:
            self.query_one("#rich_log", RichLog).write(f"\n{err_msg}")
        except Exception:
            pass


def run_tui() -> None:
    GstdTuiApp().run()


if __name__ == "__main__":
    run_tui()