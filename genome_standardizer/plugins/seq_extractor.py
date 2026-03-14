# -*- coding: utf-8 -*-

import sys
from tqdm import tqdm

try:
    from Bio.Seq import Seq
except ImportError:
    sys.exit("[FATAL] Biopython is required. Please install it: pip install biopython")


def get_sequence_slice(genome_dict, scaf, start, end, strand, missing_scafs):
    if scaf not in genome_dict:
        missing_scafs.add(scaf)
        return None

    seq_obj = genome_dict[scaf].seq[start - 1:end]
    if strand == '-':
        return seq_obj.reverse_complement()
    return seq_obj


def extract_all(genome_seqs, processed_genes, logger=None, pep_qc=False):
    cds_records = {}
    pep_records = {}
    missing_scafs    = set()
    empty_cds_count  = 0

    # QC counters (only used when --pep-qc is active)
    qc_no_atg        = 0
    qc_premature_stop = 0
    qc_truncated     = 0

    for gid, gene in tqdm(processed_genes.items(),
                          desc="         [Progress] Extracting  ",
                          unit=" gene", ncols=100, leave=False):
        scaf   = gene['scaf']
        strand = gene['strand']

        for mid, mrna in gene['mrnas'].items():
            cds_feats = [f for f in mrna['subfeats']
                         if f['type'].upper() == 'CDS']

            if not cds_feats:
                empty_cds_count += 1
                continue

            if strand == '+':
                cds_feats.sort(key=lambda x: x['start'])
            else:
                cds_feats.sort(key=lambda x: x['start'], reverse=True)

            cds_seq = Seq("")
            for cf in cds_feats:
                frag = get_sequence_slice(
                    genome_seqs, scaf, cf['start'], cf['end'],
                    strand, missing_scafs)
                if frag:
                    cds_seq += frag

            if len(cds_seq) == 0:
                continue

            cds_records[mid] = str(cds_seq)
            rem = len(cds_seq) % 3

            if pep_qc:
                # --- Translation Quality Control --------------------------
                # Track truncated CDS (not divisible by 3)
                if rem > 0:
                    qc_truncated += 1

                # Check start codon
                if not str(cds_seq[:3]).upper().startswith('ATG'):
                    qc_no_atg += 1

                # Translate full sequence to detect internal stop codons
                work_seq = cds_seq[:-rem] if rem > 0 else cds_seq
                try:
                    full_pep = str(work_seq.translate())
                    # Internal stop = stop codon before the final position
                    internal_stops = full_pep[:-1].count('*')
                    if internal_stops > 0:
                        qc_premature_stop += 1
                    # Strip trailing stop for output
                    pep_records[mid] = full_pep.rstrip('*')
                except Exception:
                    pass
            else:
                # --- Original translation (unchanged) ---------------------
                try:
                    if rem > 0:
                        pep_seq = cds_seq[:-rem].translate(to_stop=True)
                    else:
                        pep_seq = cds_seq.translate(to_stop=True)
                    pep_records[mid] = str(pep_seq)
                except Exception:
                    pass

    if missing_scafs and logger:
        sample_scafs = ", ".join(list(missing_scafs)[:5])
        logger.warning(
            f"         [Diagnostic] {len(missing_scafs)} scaffold IDs missing "
            f"in FASTA (e.g., {sample_scafs}...).")

    if empty_cds_count > 0 and logger:
        logger.warning(
            f"         [Diagnostic] {empty_cds_count} transcripts lack CDS coordinates.")

    # Report QC results
    if pep_qc and logger:
        total = len(cds_records)
        logger.info(f"\n         [PEP-QC] Translation quality report ({total} CDS total):")
        logger.info(f"         [PEP-QC]   Truncated CDS (len % 3 != 0) : {qc_truncated}")
        logger.info(f"         [PEP-QC]   Non-ATG start codon          : {qc_no_atg}")
        logger.info(f"         [PEP-QC]   Internal (premature) stop    : {qc_premature_stop}")
        if any([qc_truncated, qc_no_atg, qc_premature_stop]):
            logger.info(
                f"         [PEP-QC]   Note: All sequences are retained in output. "
                "Review flagged transcripts before downstream analysis.")

    return cds_records, pep_records
