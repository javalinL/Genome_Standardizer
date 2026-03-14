# -*- coding: utf-8 -*-

import gzip
import re
from collections import defaultdict

try:
    from Bio import SeqIO
except ImportError:
    raise ImportError("Biopython is required.")


# ---------------------------------------------------------------------------
# Lightweight proxy used by LazyFastaDict.items() — carries description only,
# no sequence data loaded.  Fully duck-type compatible with SeqRecord for the
# fuzzy-match loop (only .description and .id are accessed there).
# ---------------------------------------------------------------------------
class _RecordProxy:
    __slots__ = ('description', 'id')

    def __init__(self, description, alias):
        self.description = description
        self.id = alias


# ---------------------------------------------------------------------------
# Memory-efficient FASTA container.
# Regular mode  : genome_seqs is a plain dict {alias -> SeqRecord}.
# --low-mem mode : genome_seqs is a LazyFastaDict — stores only byte-offset
#                  indexes (SeqIO.index) + descriptions; sequences are read
#                  from disk on demand one scaffold at a time.
#                  For .gz inputs it falls back to full in-memory loading
#                  (with a warning) because gzip does not support random seek.
# ---------------------------------------------------------------------------
class LazyFastaDict:
    def __init__(self):
        self._entry   = {}   # alias -> (fasta_path_or_None, real_id_in_file)
        self._desc    = {}   # alias -> description string  (header-scan only)
        self._index   = {}   # fasta_path -> SeqIO.index handle
        self._gz_recs = {}   # real_id    -> SeqRecord  (gz fallback + fuzzy)

    def set_index(self, fpath, idx):
        self._index[fpath] = idx

    def add_full_entry(self, alias, fpath, real_id, desc=''):
        """Register a sequence ID (or alias) backed by an on-disk index."""
        self._entry[alias] = (fpath, real_id)
        self._desc[alias]  = desc

    def add_alias(self, alias, source_alias):
        """Point alias at the same entry as source_alias (no seq loaded)."""
        if source_alias in self._entry:
            self._entry[alias] = self._entry[source_alias]
            self._desc[alias]  = self._desc.get(source_alias, '')

    def __contains__(self, key):
        return key in self._entry

    def __getitem__(self, key):
        fpath, real_id = self._entry[key]
        if fpath is None or real_id in self._gz_recs:
            return self._gz_recs[real_id]
        return self._index[fpath][real_id]   # reads seq from disk on demand

    def __setitem__(self, key, value):
        """Support genome_seqs[alias] = record/proxy from the fuzzy loop."""
        if isinstance(value, _RecordProxy):
            self.add_alias(key, value.id)
        elif hasattr(value, 'seq'):            # real SeqRecord (gz fallback)
            self._gz_recs[value.id] = value
            self._entry[key] = (None, value.id)
            self._desc[key]  = getattr(value, 'description', '')

    def keys(self):
        return self._entry.keys()

    def items(self):
        """Yields (alias, _RecordProxy) — descriptions only, NO seq loading.
        Deduplicates so each underlying sequence appears once."""
        seen = set()
        for alias, (fpath, real_id) in self._entry.items():
            uid = (fpath, real_id)
            if uid not in seen:
                seen.add(uid)
                yield alias, _RecordProxy(self._desc.get(alias, ''), alias)

    def merge_from(self, other, prefix=''):
        """Absorb another LazyFastaDict, optionally prefixing all aliases."""
        for alias, (fpath, real_id) in other._entry.items():
            self._entry[f"{prefix}{alias}"] = (fpath, real_id)
            self._desc[f"{prefix}{alias}"]  = other._desc.get(alias, '')
        for fpath, idx in other._index.items():
            self._index[fpath] = idx
        for real_id, record in other._gz_recs.items():
            self._gz_recs[real_id] = record


# ---------------------------------------------------------------------------
# Chromosome name normaliser
# ---------------------------------------------------------------------------
def standardize_chr_name(name):
    name = name.strip()
    m = re.search(r'(?i)^.*?(?:chr|chromosome|lg|linkage\s*group)[_.\-\s]*([0-9]+[a-z]*)$', name)
    if m:
        num = m.group(1)
        if num.isdigit():
            num = str(int(num))
        return f"Chr{num}"
    return name


# ---------------------------------------------------------------------------
# Core per-file parser
# ---------------------------------------------------------------------------
def _parse_single_pair(gff_path, fasta_path, no_repair=False, low_mem=False):
    seq_format = "fasta"
    lower_seq_path = fasta_path.lower()
    if any(lower_seq_path.endswith(ext) for ext in
           ['.gb', '.gb.gz', '.gbk', '.gbk.gz', '.gbff', '.gbff.gz']):
        seq_format = "genbank"

    # ------------------------------------------------------------------ #
    #  SEQUENCE LOADING                                                    #
    # ------------------------------------------------------------------ #
    use_lazy = low_mem and not fasta_path.endswith('.gz') and seq_format == 'fasta'

    if use_lazy:
        # --- Low-memory path: header scan + SeqIO.index -------------------
        genome_seqs = LazyFastaDict()
        header_map  = {}   # clean_id -> (raw_id_in_file, description)

        with open(fasta_path, 'rt') as fh:
            for line in fh:
                if line.startswith('>'):
                    parts  = line[1:].strip().split(maxsplit=1)
                    raw_id = parts[0]
                    desc   = parts[1] if len(parts) > 1 else ''
                    clean_id = raw_id.split('|')[-1].strip()
                    header_map[clean_id] = (raw_id, desc)

        idx = SeqIO.index(fasta_path, 'fasta')
        genome_seqs.set_index(fasta_path, idx)

        for clean_id, (raw_id, desc) in header_map.items():
            genome_seqs.add_full_entry(clean_id, fasta_path, raw_id, desc)

            std_id = standardize_chr_name(clean_id)
            if std_id != clean_id:
                genome_seqs.add_alias(std_id, clean_id)

            chr_match = re.search(
                r'(?:chromosome|chr|linkage\s*group|LG)\s+([0-9a-zA-Z]+)',
                desc, re.IGNORECASE)
            if chr_match:
                chr_num = chr_match.group(1)
                if chr_num.isdigit():
                    chr_num = str(int(chr_num))
                genome_seqs.add_alias(f"Chr{chr_num}", clean_id)
                genome_seqs.add_alias(f"chr{chr_num}", clean_id)

            base_id   = re.sub(r'\.[0-9]+$', '', clean_id)
            num_match = re.search(r'([0-9]+)$', base_id)
            if num_match:
                full_num = num_match.group(1)
                genome_seqs.add_alias(full_num, clean_id)
                if full_num.isdigit():
                    genome_seqs.add_alias(str(int(full_num)), clean_id)

    else:
        # --- Original path: full in-memory load ---------------------------
        if low_mem and fasta_path.endswith('.gz'):
            import sys
            print(f"         [Low-Mem] Warning: {fasta_path} is gzipped — "
                  "falling back to full in-memory load. Decompress to use --low-mem.",
                  file=sys.stderr)

        handle = gzip.open(fasta_path, "rt") if fasta_path.endswith('.gz') \
            else open(fasta_path, "rt")
        genome_seqs = {}

        for record in SeqIO.parse(handle, seq_format):
            clean_id = record.id.split('|')[-1].strip()
            genome_seqs[clean_id] = record

            std_id = standardize_chr_name(clean_id)
            if std_id != clean_id:
                genome_seqs[std_id] = record

            chr_match = re.search(
                r'(?:chromosome|chr|linkage\s*group|LG)\s+([0-9a-zA-Z]+)',
                record.description, re.IGNORECASE)
            if chr_match:
                chr_num = chr_match.group(1)
                if chr_num.isdigit():
                    chr_num = str(int(chr_num))
                genome_seqs[f"Chr{chr_num}"] = record
                genome_seqs[f"chr{chr_num}"] = record

            base_id   = re.sub(r'\.[0-9]+$', '', clean_id)
            num_match = re.search(r'([0-9]+)$', base_id)
            if num_match:
                full_num = num_match.group(1)
                if full_num not in genome_seqs:
                    genome_seqs[full_num] = record
                if full_num.isdigit() and str(int(full_num)) not in genome_seqs:
                    genome_seqs[str(int(full_num))] = record

        handle.close()

    # ------------------------------------------------------------------ #
    #  GFF / GTF PARSING  (unchanged logic, works with both container types)
    # ------------------------------------------------------------------ #
    genes           = {}
    mrna_parent_map = {}
    orphan_subfeats = defaultdict(list)
    orphan_mrnas    = defaultdict(list)
    scaf_map_cache  = {}

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
                        # Fuzzy description search — works transparently for
                        # both plain dict (SeqRecord.description) and
                        # LazyFastaDict (_RecordProxy.description, no seq load)
                        for seq_id, record in genome_seqs.items():
                            if re.search(r'\b' + re.escape(raw_scaf) + r'\b',
                                         record.description) or \
                               re.search(r'\b' + re.escape(scaf) + r'\b',
                                         record.description):
                                genome_seqs[raw_scaf] = record
                                genome_seqs[scaf]     = record
                                scaf = raw_scaf
                                found_match = True
                                break
                scaf_map_cache[raw_scaf] = scaf

            start, end = int(parts[3]), int(parts[4])
            score, strand, phase, attr_str = parts[5], parts[6], parts[7], parts[8]

            attrs  = {}
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
                        'scaf': scaf, 'start': start, 'end': end,
                        'strand': strand, 'source': parts[1], 'attrs': attrs,
                        'mrnas': defaultdict(list)
                    }
                else:
                    genes[gid]['start'] = min(genes[gid]['start'], start)
                    genes[gid]['end']   = max(genes[gid]['end'],   end)
                if mid:
                    if mid not in genes[gid]['mrnas']:
                        genes[gid]['mrnas'][mid] = []
                        mrna_parent_map[mid] = gid
                    if ftype in ['CDS', 'exon', 'five_prime_UTR',
                                 'three_prime_UTR', 'UTR']:
                        genes[gid]['mrnas'][mid].append(feat)
            else:
                if 'ID' not in attrs and 'Parent' not in attrs: continue
                if ftype == 'gene':
                    gid = attrs['ID']
                    genes[gid] = {
                        'scaf': scaf, 'start': start, 'end': end,
                        'strand': strand, 'source': parts[1], 'attrs': attrs,
                        'mrnas': defaultdict(list)
                    }
                elif ftype in ['mRNA', 'transcript']:
                    mid  = attrs['ID']
                    pids = attrs.get('Parent', '').split(',')
                    for pid in pids:
                        if not pid: continue
                        if pid in genes:
                            genes[pid]['mrnas'][mid] = []
                            mrna_parent_map[mid] = pid
                        else:
                            orphan_mrnas[pid].append((mid, feat))
                elif ftype in ['CDS', 'exon', 'five_prime_UTR',
                               'three_prime_UTR', 'UTR']:
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
            gid       = pid
            dummy_mid = f"{gid}_transcript_1"
            if dummy_mid not in genes[gid]['mrnas']:
                genes[gid]['mrnas'][dummy_mid] = []
                mrna_parent_map[dummy_mid] = gid
            genes[gid]['mrnas'][dummy_mid].extend(feats)

    # ------------------------------------------------------------------ #
    #  AUTO-HEAL  (skipped when --no-repair is set)                        #
    # ------------------------------------------------------------------ #
    repair_counts = {'cds_from_exon': 0, 'exon_from_cds': 0}

    if not no_repair:
        for gid, gene_data in genes.items():
            for mid, subfeats in gene_data['mrnas'].items():
                has_cds  = any(f['type'].upper() == 'CDS'  for f in subfeats)
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


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------
def parse_inputs(gff_list, fasta_list, custom_prefixes=None, logger=None,
                 no_repair=False, low_mem=False):
    global_genome_seqs = LazyFastaDict() if low_mem else {}
    global_genes       = {}

    for idx, (gff_path, fasta_path) in enumerate(zip(gff_list, fasta_list)):
        sub_prefix = custom_prefixes[idx] if custom_prefixes else ""
        file_name  = gff_path.split('/')[-1]

        if logger:
            logger.info(f"         --> Parsing chunk {idx + 1}: {file_name} ...")

        local_seqs, local_genes, repair_counts = _parse_single_pair(
            gff_path, fasta_path, no_repair=no_repair, low_mem=low_mem)

        if logger:
            if repair_counts['cds_from_exon'] > 0:
                logger.info(
                    f"         [Auto-Heal] Reconstructed CDS for "
                    f"{repair_counts['cds_from_exon']} mRNAs based on Exon coordinates.")
            if repair_counts['exon_from_cds'] > 0:
                logger.info(
                    f"         [Auto-Heal] Reconstructed Exon for "
                    f"{repair_counts['exon_from_cds']} mRNAs based on CDS coordinates.")

        # Merge sequences
        if low_mem:
            global_genome_seqs.merge_from(local_seqs, prefix=sub_prefix)
            if not custom_prefixes:
                for seq_key in local_seqs.keys():
                    new_key = f"{sub_prefix}{seq_key}"
                    if new_key in global_genome_seqs:
                        fpath_a, rid_a = global_genome_seqs._entry[new_key]
                        fpath_b, rid_b = local_seqs._entry.get(seq_key, (None, None))
                        if (fpath_a, rid_a) != (fpath_b, rid_b):
                            raise ValueError(
                                f"\n[Conflict Radar] Fatal sequence conflict!\n"
                                f"Scaffold '{new_key}' exists in multiple chunks.\n"
                                f"Fix: Use --add-prefix (e.g., SubC_,SubD_) to isolate subgenomes.")
        else:
            for seq_key, record in local_seqs.items():
                new_seq_key = f"{sub_prefix}{seq_key}"
                if not custom_prefixes and new_seq_key in global_genome_seqs:
                    raise ValueError(
                        f"\n[Conflict Radar] Fatal sequence conflict!\n"
                        f"Scaffold '{new_seq_key}' exists in multiple chunks.\n"
                        f"Fix: Use --add-prefix (e.g., SubC_,SubD_) to isolate subgenomes.")
                global_genome_seqs[new_seq_key] = record

        # Merge genes
        for gid, gene_data in local_genes.items():
            new_gid = f"{sub_prefix}{gid}"
            if not custom_prefixes and new_gid in global_genes:
                raise ValueError(
                    f"\n[Conflict Radar] Fatal gene conflict!\n"
                    f"Gene '{new_gid}' exists in multiple chunks.\n"
                    f"Fix: Use --add-prefix (e.g., SubC_,SubD_) to isolate subgenomes.")

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
