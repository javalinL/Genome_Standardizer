# -*- coding: utf-8 -*-

import sys

try:
    from Bio.Seq import Seq
except ImportError:
    sys.exit("[FATAL] Biopython is required. Please install it: pip install biopython")

def get_sequence_slice(genome_dict, scaf, start, end, strand, missing_scafs, logger):
    if scaf not in genome_dict:
        if scaf not in missing_scafs:
            if logger:
                logger.warning(f"\n[Warning] Sequence mismatch! Scaffold '{scaf}' found in GFF but missing in FASTA.")
            missing_scafs.add(scaf)
        return None

    seq_obj = genome_dict[scaf].seq[start - 1:end]
    if strand == '-':
        return seq_obj.reverse_complement()
    return seq_obj

def extract_all(genome_seqs, processed_genes, logger=None):
    cds_records = {}
    pep_records = {}
    missing_scafs = set()
    empty_cds_count = 0

    for gid, gene in processed_genes.items():
        scaf = gene['scaf']
        strand = gene['strand']

        for mid, mrna in gene['mrnas'].items():
            cds_feats = [f for f in mrna['subfeats'] if f['type'].upper() == 'CDS']

            if not cds_feats:
                empty_cds_count += 1
                continue

            if strand == '+':
                cds_feats.sort(key=lambda x: x['start'])
            else:
                cds_feats.sort(key=lambda x: x['start'], reverse=True)

            cds_seq = Seq("")
            for cf in cds_feats:
                frag = get_sequence_slice(genome_seqs, scaf, cf['start'], cf['end'], strand, missing_scafs, logger)
                if frag:
                    cds_seq += frag

            if len(cds_seq) > 0:
                cds_records[mid] = str(cds_seq)
                rem = len(cds_seq) % 3
                try:
                    if rem > 0:
                        pep_seq = cds_seq[:-rem].translate(to_stop=True)
                    else:
                        pep_seq = cds_seq.translate(to_stop=True)
                    pep_records[mid] = str(pep_seq)
                except Exception:
                    pass

    if missing_scafs and logger:
        logger.warning(f"         [Diagnostic] {len(missing_scafs)} scaffold IDs missing in FASTA.")
    if empty_cds_count > 0 and logger:
        logger.warning(f"         [Diagnostic] {empty_cds_count} transcripts lack CDS coordinates.")

    return cds_records, pep_records