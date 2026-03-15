# Genome Standardizer (`gstd`) v4.0.0

[![Python 3.7+](https://img.shields.io/badge/python-3.7%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Textual TUI](https://img.shields.io/badge/TUI-Textual%208.x-purple.svg)](https://github.com/Textualize/textual)

A robust, memory-efficient pipeline for standardizing plant genome annotations, resolving structural inconsistencies, and extracting high-quality coding sequences for comparative genomics. v4.0.0 introduces an optional, fully interactive **Terminal User Interface (TUI)** that operates natively over SSH without X11 forwarding, making it suitable for remote HPC environments.

---

## Overview

Working with publicly available genome assemblies frequently involves heterogeneous GFF3/GTF files with broken feature hierarchies, mismatched sequence identifiers between FASTA headers and coordinate lines, and non-standard attribute formatting. `gstd` acts as a universal annotation adapter: it takes raw FASTA and GFF3/GTF files as input and produces a coordinated, strictly structured output set ready for downstream phylogenetic inference, synteny analysis, and evolutionary rate (Ka/Ks) calculation.

---

## Installation

Python 3.7+ is required. Two installation modes are available to keep dependencies minimal when the TUI is not needed.

**Full TUI mode** (recommended — installs the Textual interactive frontend):

```bash
pip install "genome-standardizer[tui] @ git+https://github.com/javalinL/Genome_Standardizer.git"
```

**Minimal CLI mode** (core pipeline only, no additional UI dependencies):

```bash
pip install git+https://github.com/javalinL/Genome_Standardizer.git
```

---

## Quick Start

### Interactive TUI

Launch the graphical terminal interface. Use `Ctrl+L` to switch between English and Chinese; press `F1` to open the integrated reference manual.

```bash
gstd-tui
```

### Headless CLI

The `gstd` command is available after installation regardless of mode.

**Standard execution:**

```bash
gstd annotation.gff3 genome.fasta Target_Prefix --save-log
```

**Longest-transcript mode** — recommended for phylogenetic tree construction and Ka/Ks calculations:

```bash
gstd annotation.gff3 genome.fasta Target_Prefix --longest --save-log
```

**Polyploid subgenome fusion** — assigns isolated ID prefixes to prevent fatal collision during merge:

```bash
gstd subA.gff,subB.gff subA.fa,subB.fa Hybrid_Species --add-prefix SubA_,SubB_ --longest
```

**Large genome mode** — recommended for assemblies ≥ 10 GB (e.g., hexaploid wheat, sugarcane); requires uncompressed FASTA:

```bash
gstd annotation.gff3 genome.fasta Target_Prefix --low-mem --save-log
```

**Strict public-database input with full QC and provenance tracking:**

```bash
gstd annotation.gff3 genome.fasta Target_Prefix --no-repair --pep-qc --keep-source --bed6 --save-log
```

> Use `--no-repair` when inputs contain non-coding gene models (lncRNA, tRNA, pseudogenes) annotated with exon-only features; the default Auto-Heal mode would otherwise infer spurious CDS for such entries. Use `--pep-qc` to assess protein sequence integrity prior to alignment or clustering.

---

## Key Features

| Feature | Description |
|---|---|
| **Interactive TUI** | Mouse-supported bilingual (EN/ZH) terminal interface with real-time pipeline progress. No X11 required; fully functional over SSH. |
| **HPC Job Script Generator** | Built-in SLURM and PBS/LSF script generator. Cluster resources and scheduler directives are configured visually and exported as a ready-to-submit shell script. |
| **Smart Alias Sniffer** | Automatically detects and resolves mismatches between FASTA sequence identifiers and GFF3 coordinate-line sequence names. |
| **Annotation Auto-Heal** | Reconstructs broken GFF hierarchies: infers missing CDS from exon coordinates, generates parent `gene` features for orphan transcripts. Disable with `--no-repair` for non-coding gene sets. |
| **Polyploid Collision Radar** | Detects sequence and gene ID collisions before subgenome fusion, preventing silent data overwriting in allopolyploid assemblies. |
| **Low-Memory Index Mode** | On-disk byte-offset FASTA indexing via `SeqIO.index()`. Sequences are loaded one scaffold at a time, eliminating OOM failures on assemblies such as hexaploid wheat (≥ 15 GB) or sugarcane. Requires uncompressed FASTA input. |
| **PEP Quality Control** | Reports counts of non-ATG start codons, internal premature stop codons, and truncated CDS at run completion. All sequences are retained in output; results are diagnostic only. |
| **Provenance Preservation** | Optionally retains the original annotation source tag (e.g., `augustus`, `maker`, `SNAP`) in GFF3 column 2 (`--keep-source`). |

---

## CLI Options Reference

| Flag | Default | Description |
|---|---|---|
| `--step N` | `10` | Numeric increment between consecutive standardized gene IDs |
| `--add-prefix X,Y,...` | — | Assign isolated ID prefixes to each subgenome for polyploid fusion |
| `--longest` | off | Retain only the longest transcript (representative isoform) per locus |
| `--keep` | off | Skip gzip compression of input files after processing |
| `--save-log` | off | Write a timestamped `<Prefix>_run.log` file to the working directory |
| `--no-repair` | off | Disable Auto-Heal CDS/exon inference |
| `--pep-qc` | off | Report translation quality statistics at run completion |
| `--bed6` | off | Output strict 6-column BED format (compatible with `bedtools`, IGV) |
| `--keep-source` | off | Preserve original annotation source field in GFF3 column 2 |
| `--low-mem` | off | Use on-disk FASTA indexing instead of loading sequences into RAM |

---

## Output Files

Five standardized files are produced in the working directory for each processed species:

| File | Description |
|---|---|
| `ALL.standard.gff3` | Structurally correct GFF3 with strictly sorted gene → mRNA → CDS/exon hierarchies |
| `<Prefix>.bed` | Gene coordinate file. Extended 7-column format by default (column 7 = original ID); use `--bed6` for standard 6-column BED |
| `<Prefix>.cds` | FASTA file of coding sequences, one entry per representative transcript |
| `<Prefix>.pep` | Translated protein sequences, ready for alignment or orthology inference |
| `<Prefix>.id_list` | Tab-delimited mapping table from original annotation IDs to standardized IDs |

---

## Acknowledgements

The v4.0.0 release was developed using a structured Dual-Agent AI workflow under the direction of the lead developer:

- **L. javalin** (Lead Developer & Project Owner): Designed the core algorithmic logic, defined the bioinformatic standards, and orchestrated the multi-agent collaboration.
- **Claude Sonnet 4.6** (Anthropic): Engineered the asynchronous pipeline wrapper, thread-safe message queue, i18n architecture, HPC Job Script Generator, and integration test suite.
- **Gemini 2.5 Pro** (Google DeepMind): Designed and implemented the complete Textual UI layout (`app.py`), dual-theme stylesheet (`app.tcss`), and bioinformatically refined English/Chinese localizations.

---

## License

MIT License. Created by L. javalin.
