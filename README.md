# Genome Standardizer (gstd)

A highly robust, memory-efficient pipeline for standardizing complex plant genomes, resolving annotation inconsistencies, and extracting high-quality sequences for comparative genomics.

## Overview

Working with publicly available genome assemblies often involves dealing with highly fragmented scaffolds, missing hierarchical annotations (e.g., missing `gene` or `exon` features), and conflicting sequence IDs. `gstd` acts as a universal adapter, taking raw FASTA and GFF3/GTF files and producing clean, standardized, and strictly coordinated outputs ready for downstream phylogenetic and synteny analyses.

## Key Features

* **Smart Alias Sniffer**: Automatically detects and heals mismatching sequence IDs between FASTA headers and GFF coordinate lines.
* **Annotation Auto-Healing**: Reconstructs broken GFF hierarchies (e.g., inferring missing CDS from exons, or generating parent `gene` features for orphan transcripts).
* **Polyploid Fusion & Collision Radar**: Safely merges subgenomes of complex polyploids. The built-in radar prevents fatal data overwriting by detecting sequence or gene ID collisions before they corrupt the database.
* **Longest Transcript Selection**: Optional filtering to retain only the most representative isoform per locus, reducing redundancy for gene family clustering.
* **Extreme Memory Efficiency**: Stream-based sequence processing prevents Out-Of-Memory (OOM) kills, even when handling assemblies with millions of unanchored contigs.

## Installation

Ensure you have Python 3.7+ and Biopython installed. You can install `gstd` globally via pip directly from GitHub:

pip install git+https://github.com/javalinL/Genome_Standardizer.git

## Quick Start
The installation automatically binds the tool to the gstd command globally.

1. Standard Execution (All Isoforms):

gstd annotation.gff3 genome.fasta Target_Prefix --save-log

2. Retain Longest Transcript Only:Highly recommended for phylogenetic tree construction and $Ka/Ks$ calculations.

gstd annotation.gff3 genome.fasta Target_Prefix --longest --save-log

3. Polyploid Subgenome Fusion:
Merge subgenomes safely by assigning isolated prefixes to prevent ID collisions:

gstd subA.gff,subB.gff subA.fa,subB.fa Hybrid_Species --add-prefix SubA_,SubB_ --save-log

## Output Files
The pipeline generates five standardized files in the working directory:

1.ALL.standard.gff3: Structurally perfect annotation file with strictly sorted hierarchies.

2.<Prefix>.bed: Lightweight coordinate file for fast interval operations (e.g., bedtools).

3.<Prefix>.cds: Clean nucleotide sequences for coding regions.

4.<Prefix>.pep: Translated protein sequences ready for alignment.

5.<Prefix>.id_list: A mapping dictionary bridging original erratic IDs to the new, clean standardized IDs.

## Acknowledgements
This tool was collaboratively developed with the assistance of Gemini (Google's AI model), serving as an AI pair programmer. Gemini contributed to the core logic refactoring, robust error-handling mechanisms (including the Collision Radar and Smart Alias Sniffer), and the English localization of the codebase and documentation.

## License
MIT License. Created by L. javalin.

