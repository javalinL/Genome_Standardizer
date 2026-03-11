# -*- coding: utf-8 -*-

import gzip
import re
from collections import defaultdict

try:
    from Bio import SeqIO
except ImportError:
    raise ImportError("Biopython is required.")


def standardize_chr_name(name):
    name = name.strip()
    m = re.search(r'(?i)^.*?(?:chr|chromosome|lg|linkage\s*group)[_.\-\s]*([0-9]+[a-z]*)$', name)
    if m:
        num = m.group(1)
        if num.isdigit():
            num = str(int(num))
        return f"Chr{num}"
    return name


def _parse_single_pair(gff_path, fasta_path):
    seq_format = "fasta"
    lower_seq_path = fasta_path.lower()
    if any(lower_seq_path.endswith(ext) for ext in ['.gb', '.gb.gz', '.gbk', '.gbk.gz', '.gbff', '.gbff.gz']):
        seq_format = "genbank"

    handle = gzip.open(fasta_path, "rt") if fasta_path.endswith('.gz') else open(fasta_path, "rt")
    genome_seqs = {}

    for record in SeqIO.parse(handle, seq_format):
        clean_id = record.id.split('|')[-1].strip()
        genome_seqs[clean_id] = record

        std_id = standardize_chr_name(clean_id)
        if std_id != clean_id:
            genome_seqs[std_id] = record

        chr_match = re.search(r'(?:chromosome|chr|linkage\s*group|LG)\s+([0-9a-zA-Z]+)', record.description,
                              re.IGNORECASE)
        if chr_match:
            chr_num = chr_match.group(1)
            if chr_num.isdigit():
                chr_num = str(int(chr_num))
            genome_seqs[f"Chr{chr_num}"] = record
            genome_seqs[f"chr{chr_num}"] = record

        base_id = re.sub(r'\.[0-9]+$', '', clean_id)
        num_match = re.search(r'([0-9]+)$', base_id)
        if num_match:
            full_num = num_match.group(1)
            if full_num not in genome_seqs:
                genome_seqs[full_num] = record
            if full_num.isdigit() and str(int(full_num)) not in genome_seqs:
                genome_seqs[str(int(full_num))] = record

    handle.close()

    genes = {}
    mrna_parent_map = {}
    orphan_subfeats = defaultdict(list)
    orphan_mrnas = defaultdict(list)
    scaf_map_cache = {}

    open_func = gzip.open if gff_path.endswith('.gz') else open

    with open_func(gff_path, 'rt') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) != 9: continue

            raw_scaf = parts[0]
            if raw_scaf in scaf_map_cache:
                scaf = scaf_map_cache[raw_scaf]
            else:
                scaf = standardize_chr_name(raw_scaf)
                if scaf not in genome_seqs and raw_scaf not in genome_seqs:
                    found_match = False
                    for target in [scaf, raw_scaf]:
                        num_match = re.search(r'([0-9]+)$', target)
                        if num_match:
                            full_num = num_match.group(1)
                            if full_num in genome_seqs:
                                scaf = full_num
                                found_match = True
                                break
                            elif str(int(full_num)) in genome_seqs:
                                scaf = str(int(full_num))
                                found_match = True
                                break
                    if not found_match:
                        for seq_id, record in genome_seqs.items():
                            if re.search(r'\b' + re.escape(raw_scaf) + r'\b', record.description) or \
                                    re.search(r'\b' + re.escape(scaf) + r'\b', record.description):
                                genome_seqs[raw_scaf] = record
                                genome_seqs[scaf] = record
                                scaf = raw_scaf
                                found_match = True
                                break
                scaf_map_cache[raw_scaf] = scaf

            start, end = int(parts[3]), int(parts[4])
            score, strand, phase, attr_str = parts[5], parts[6], parts[7], parts[8]

            attrs = {}
            is_gtf = False
            if '=' in attr_str:
                for item in attr_str.strip().split(';'):
                    if not item: continue
                    if '=' in item:
                        k, v = item.split('=', 1)
                        attrs[k.strip()] = v.strip()
            else:
                is_gtf = True
                for m in re.finditer(r'([a-zA-Z0-9_.-]+)\s+"([^"]+)"', attr_str):
                    attrs[m.group(1)] = m.group(2)

            feat = {
                'scaf': scaf, 'source': parts[1], 'type': parts[2],
                'start': start, 'end': end, 'score': score,
                'strand': strand, 'phase': phase, 'attrs': attrs
            }
            ftype = parts[2]

            if is_gtf:
                gid = attrs.get('gene_id')
                mid = attrs.get('transcript_id')
                if not gid: continue
                if gid not in genes:
                    genes[gid] = {
                        'scaf': scaf, 'start': start, 'end': end, 'strand': strand,
                        'source': parts[1], 'attrs': attrs,
                        'mrnas': defaultdict(list)
                    }
                else:
                    genes[gid]['start'] = min(genes[gid]['start'], start)
                    genes[gid]['end'] = max(genes[gid]['end'], end)
                if mid:
                    if mid not in genes[gid]['mrnas']:
                        genes[gid]['mrnas'][mid] = []
                        mrna_parent_map[mid] = gid
                    if ftype in ['CDS', 'exon', 'five_prime_UTR', 'three_prime_UTR', 'UTR']:
                        genes[gid]['mrnas'][mid].append(feat)
            else:
                if 'ID' not in attrs and 'Parent' not in attrs: continue
                if ftype == 'gene':
                    gid = attrs['ID']
                    genes[gid] = {
                        'scaf': scaf, 'start': start, 'end': end, 'strand': strand,
                        'source': parts[1], 'attrs': attrs,
                        'mrnas': defaultdict(list)
                    }
                elif ftype in ['mRNA', 'transcript']:
                    mid = attrs['ID']
                    pids = attrs.get('Parent', '').split(',')
                    for pid in pids:
                        if not pid: continue
                        if pid in genes:
                            genes[pid]['mrnas'][mid] = []
                            mrna_parent_map[mid] = pid
                        else:
                            orphan_mrnas[pid].append((mid, feat))
                elif ftype in ['CDS', 'exon', 'five_prime_UTR', 'three_prime_UTR', 'UTR']:
                    pids = attrs.get('Parent', '').split(',')
                    for pid in pids:
                        if not pid: continue
                        if pid in mrna_parent_map:
                            gid = mrna_parent_map[pid]
                            genes[gid]['mrnas'][pid].append(feat)
                        else:
                            orphan_subfeats[pid].append(feat)

    for pid, mrna_list in orphan_mrnas.items():
        if pid in genes:
            for mid, mfeat in mrna_list:
                if mid not in genes[pid]['mrnas']:
                    genes[pid]['mrnas'][mid] = []
                mrna_parent_map[mid] = pid

    for pid, feats in orphan_subfeats.items():
        if pid in mrna_parent_map:
            gid = mrna_parent_map[pid]
            genes[gid]['mrnas'][pid].extend(feats)
        elif pid in genes:
            gid = pid
            dummy_mid = f"{gid}_transcript_1"
            if dummy_mid not in genes[gid]['mrnas']:
                genes[gid]['mrnas'][dummy_mid] = []
                mrna_parent_map[dummy_mid] = gid
            genes[gid]['mrnas'][dummy_mid].extend(feats)

    repair_counts = {'cds_from_exon': 0, 'exon_from_cds': 0}
    for gid, gene_data in genes.items():
        for mid, subfeats in gene_data['mrnas'].items():
            has_cds = any(f['type'].upper() == 'CDS' for f in subfeats)
            has_exon = any(f['type'].upper() == 'EXON' for f in subfeats)
            if has_exon and not has_cds:
                new_cds_feats = []
                for f in subfeats:
                    if f['type'].upper() == 'EXON':
                        new_cds = f.copy()
                        new_cds['type'] = 'CDS'
                        new_cds_feats.append(new_cds)
                subfeats.extend(new_cds_feats)
                repair_counts['cds_from_exon'] += 1
            elif has_cds and not has_exon:
                new_exon_feats = []
                for f in subfeats:
                    if f['type'].upper() == 'CDS':
                        new_exon = f.copy()
                        new_exon['type'] = 'exon'
                        new_exon_feats.append(new_exon)
                subfeats.extend(new_exon_feats)
                repair_counts['exon_from_cds'] += 1

    return genome_seqs, genes, repair_counts


def parse_inputs(gff_list, fasta_list, custom_prefixes=None, logger=None):
    global_genome_seqs = {}
    global_genes = {}

    for idx, (gff_path, fasta_path) in enumerate(zip(gff_list, fasta_list)):
        sub_prefix = custom_prefixes[idx] if custom_prefixes else ""
        file_name = gff_path.split('/')[-1]

        if logger:
            logger.info(f"         --> Parsing chunk {idx + 1}: {file_name} ...")

        local_seqs, local_genes, repair_counts = _parse_single_pair(gff_path, fasta_path)

        if logger:
            if repair_counts['cds_from_exon'] > 0:
                logger.info(
                    f"         [Auto-Heal] Reconstructed CDS for {repair_counts['cds_from_exon']} mRNAs based on Exon coordinates.")
            if repair_counts['exon_from_cds'] > 0:
                logger.info(
                    f"         [Auto-Heal] Reconstructed Exon for {repair_counts['exon_from_cds']} mRNAs based on CDS coordinates.")

        for seq_key, record in local_seqs.items():
            new_seq_key = f"{sub_prefix}{seq_key}"
            if not custom_prefixes and new_seq_key in global_genome_seqs:
                raise ValueError(
                    f"\n[Conflict Radar] Fatal sequence conflict!\n"
                    f"Scaffold '{new_seq_key}' exists in multiple chunks.\n"
                    f"Fix: Use --add-prefix (e.g., SubC_,SubD_) to isolate subgenomes."
                )
            global_genome_seqs[new_seq_key] = record

        for gid, gene_data in local_genes.items():
            new_gid = f"{sub_prefix}{gid}"
            if not custom_prefixes and new_gid in global_genes:
                raise ValueError(
                    f"\n[Conflict Radar] Fatal gene conflict!\n"
                    f"Gene '{new_gid}' exists in multiple chunks.\n"
                    f"Fix: Use --add-prefix (e.g., SubC_,SubD_) to isolate subgenomes."
                )

            gene_data['scaf'] = f"{sub_prefix}{gene_data['scaf']}"

            new_mrnas = {}
            for mid, subfeats in gene_data['mrnas'].items():
                new_mid = f"{sub_prefix}{mid}"
                for feat in subfeats:
                    feat['scaf'] = f"{sub_prefix}{feat['scaf']}"
                new_mrnas[new_mid] = subfeats
            gene_data['mrnas'] = new_mrnas

            global_genes[new_gid] = gene_data

    return global_genome_seqs, global_genes