# -*- coding: utf-8 -*-

import re
from tqdm import tqdm


def natural_sort_key(s):
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split('([0-9]+)', s)]


def _calculate_transcript_length(subfeatures):
    cds_feats = [f for f in subfeatures if f['type'].upper() == 'CDS']
    if cds_feats:
        return sum(f['end'] - f['start'] + 1 for f in cds_feats)

    exon_feats = [f for f in subfeatures if f['type'].upper() == 'EXON']
    if exon_feats:
        return sum(f['end'] - f['start'] + 1 for f in exon_feats)

    starts = [f['start'] for f in subfeatures]
    ends   = [f['end']   for f in subfeatures]
    if starts and ends:
        return max(ends) - min(starts) + 1
    return 0


def standardize_and_rename(raw_genes, prefix, step=10,
                            longest_only=False, keep_source=False):
    processed_genes = {}
    id_mapping      = []

    sorted_gids = sorted(
        raw_genes.keys(),
        key=lambda k: (natural_sort_key(raw_genes[k]['scaf']),
                       raw_genes[k]['start'])
    )

    counter = 0

    for gid in tqdm(sorted_gids,
                    desc="         [Progress] Structuring ",
                    unit=" gene", ncols=100, leave=False):
        gene  = raw_genes[gid]
        mrnas = gene.get('mrnas', {})

        if not mrnas:
            continue

        if longest_only and len(mrnas) > 1:
            longest_mid = max(mrnas.keys(),
                              key=lambda mid: _calculate_transcript_length(mrnas[mid]))
            mrnas = {longest_mid: mrnas[longest_mid]}

        counter += step
        new_gene_id = f"{prefix}_{counter:06d}"

        gene_entry = {
            'scaf':   gene['scaf'],
            'start':  gene['start'],
            'end':    gene['end'],
            'strand': gene['strand'],
            'attrs':  {'ID': new_gene_id, 'Name': new_gene_id},
            'mrnas':  {}
        }
        # Preserve original source annotation when requested
        if keep_source:
            gene_entry['source'] = gene.get('source', '.')

        processed_genes[new_gene_id] = gene_entry

        sorted_mids = sorted(mrnas.keys())

        for idx, old_mid in enumerate(sorted_mids, 1):
            new_mrna_id = f"{prefix}.{counter:06d}.{idx}"
            id_mapping.append(f"{old_mid}\t{new_mrna_id}")

            subfeatures = mrnas[old_mid]
            all_starts  = [f['start'] for f in subfeatures]
            all_ends    = [f['end']   for f in subfeatures]
            m_start = min(all_starts) if all_starts else gene['start']
            m_end   = max(all_ends)   if all_ends   else gene['end']

            new_mrna = {
                'start': m_start,
                'end':   m_end,
                'attrs': {'ID': new_mrna_id, 'Parent': new_gene_id,
                          'Name': new_mrna_id},
                'subfeats': []
            }

            subfeatures.sort(key=lambda x: x['start'])
            type_counts = {}
            for f in subfeatures:
                ftype = f['type']
                type_counts[ftype] = type_counts.get(ftype, 0) + 1
                sub_id = f"{new_mrna_id}.{ftype}{type_counts[ftype]}"

                new_feat = {
                    'type':  ftype,
                    'start': f['start'],
                    'end':   f['end'],
                    'phase': f['phase'],
                    'attrs': {'ID': sub_id, 'Parent': new_mrna_id}
                }
                # Preserve original source per sub-feature when requested
                if keep_source:
                    new_feat['source'] = f.get('source', '.')

                new_mrna['subfeats'].append(new_feat)

            processed_genes[new_gene_id]['mrnas'][new_mrna_id] = new_mrna

    return processed_genes, id_mapping
