# Genome Standardizer (gstd) v3.0.0

A highly robust, memory-efficient pipeline for standardizing complex plant genomes, resolving annotation inconsistencies, and extracting high-quality sequences for comparative genomics.

## Overview

Working with publicly available genome assemblies often involves dealing with highly fragmented scaffolds, missing hierarchical annotations (e.g., missing `gene` or `exon` features), and conflicting sequence IDs. `gstd` acts as a universal adapter, taking raw FASTA and GFF3/GTF files and producing clean, standardized, and strictly coordinated outputs ready for downstream phylogenetic and synteny analyses.

## Key Features

* **Smart Alias Sniffer**: Automatically detects and heals mismatching sequence IDs between FASTA headers and GFF coordinate lines.
* **Annotation Auto-Healing**: Reconstructs broken GFF hierarchies (e.g., inferring missing CDS from exons, or generating parent `gene` features for orphan transcripts). Can be disabled with `--no-repair` when inputs contain non-coding genes (lncRNA, tRNA, pseudogenes) to prevent false CDS generation.
* **Polyploid Fusion & Collision Radar**: Safely merges subgenomes of complex polyploids. The built-in radar prevents fatal data overwriting by detecting sequence or gene ID collisions before they corrupt the database.
* **Longest Transcript Selection**: Optional filtering to retain only the most representative isoform per locus, reducing redundancy for gene family clustering.
* **Large Genome Memory Mode (`--low-mem`)**: On-disk byte-offset indexing via `SeqIO.index()` — sequences are loaded one scaffold at a time instead of all at once. Eliminates OOM kills on assemblies like hexaploid wheat (15 GB) or sugarcane. Requires uncompressed FASTA input.
* **Translation Quality Control (`--pep-qc`)**: Reports counts of non-ATG start codons, internal (premature) stop codons, and truncated CDS at run completion. All sequences are retained in output; flags are for diagnostics only.
* **Standard BED6 Output (`--bed6`)**: Outputs strict 6-column BED format (chrom/start/end/name/score/strand) for full compatibility with `bedtools`, IGV, and genome browsers. Default retains extended 7-column format with original ID mapping in column 7.
* **Source Field Preservation (`--keep-source`)**: Retains the original annotation source (e.g., `augustus`, `maker`, `SNAP`) in GFF3 output column 2 for provenance tracking. Default outputs `.` for all features.

## Installation

Ensure you have Python 3.7+ and Biopython installed. You can install `gstd` globally via pip directly from GitHub:

```
pip install git+https://github.com/javalinL/Genome_Standardizer.git
```

## Quick Start

The installation automatically binds the tool to the `gstd` command globally.

**1. Standard Execution (All Isoforms):**
```
gstd annotation.gff3 genome.fasta Target_Prefix --save-log
```

**2. Retain Longest Transcript Only:**
Highly recommended for phylogenetic tree construction and Ka/Ks calculations.
```
gstd annotation.gff3 genome.fasta Target_Prefix --longest --save-log
```

**3. Polyploid Subgenome Fusion:**
Merge subgenomes safely by assigning isolated prefixes to prevent ID collisions:
```
gstd subA.gff,subB.gff subA.fa,subB.fa Hybrid_Species --add-prefix SubA_,SubB_ --save-log
```

**4. Large Genome (e.g., wheat, sugarcane) — prevent OOM:**
```
gstd annotation.gff3 genome.fasta Target_Prefix --low-mem --save-log
```

**5. Strict public genome input (with QC and provenance):**
```
gstd annotation.gff3 genome.fasta Target_Prefix --no-repair --pep-qc --keep-source --bed6 --save-log
```
> Use `--no-repair` when the input genome contains non-coding RNA genes or pseudogenes annotated with exon-only features; without it, those entries will have CDS incorrectly inferred. Use `--pep-qc` to assess protein sequence integrity before downstream analyses.

## Options Reference

| Flag | Description |
|---|---|
| `--step N` | Step size for gene numbering (default: 10) |
| `--add-prefix X,Y` | Assign isolated prefixes to subgenomes for polyploid fusion |
| `--longest` | Keep only the longest transcript per gene locus |
| `--keep` | Skip compression of original input files |
| `--save-log` | Write a detailed timestamped log file |
| `--no-repair` | Disable Auto-Heal CDS/exon inference |
| `--pep-qc` | Report translation quality statistics |
| `--bed6` | Output standard 6-column BED format |
| `--keep-source` | Preserve original source field in GFF3 output |
| `--low-mem` | Use on-disk indexing for memory-efficient FASTA loading |

## Output Files

The pipeline generates five standardized files in the working directory:

1. `ALL.standard.gff3`: Structurally perfect annotation file with strictly sorted hierarchies.
2. `<Prefix>.bed`: Coordinate file for interval operations. Extended 7-col format by default (col 7 = original ID); use `--bed6` for standard 6-col BED.
3. `<Prefix>.cds`: Clean nucleotide sequences for coding regions.
4. `<Prefix>.pep`: Translated protein sequences ready for alignment.
5. `<Prefix>.id_list`: A mapping dictionary bridging original erratic IDs to the new, clean standardized IDs.

## Acknowledgements

This tool was collaboratively developed with the assistance of multiple AI systems serving as pair programmers:

- **Gemini** (Google): contributed to core logic refactoring, robust error-handling mechanisms (Collision Radar, Smart Alias Sniffer), and English localization of the codebase.
- **Claude Sonnet 4.6** (Anthropic): contributed the v3.0.0 major feature update — including the `LazyFastaDict` memory-efficient FASTA architecture (`--low-mem`), translation quality control (`--pep-qc`), Auto-Heal control (`--no-repair`), standard BED6 output (`--bed6`), and source field preservation (`--keep-source`).

## License

MIT License. Created by L. javalin.
