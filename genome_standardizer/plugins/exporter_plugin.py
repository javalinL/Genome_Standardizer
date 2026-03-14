# -*- coding: utf-8 -*-

import os
from tqdm import tqdm


def export_all(processed_genes, cds_records, pep_records, id_mapping,
               prefix, work_dir, logger=None,
               bed6=False, keep_source=False):

    out_gff = os.path.join(work_dir, "ALL.standard.gff3")
    out_bed = os.path.join(work_dir, f"{prefix}.bed")
    out_cds = os.path.join(work_dir, f"{prefix}.cds")
    out_pep = os.path.join(work_dir, f"{prefix}.pep")
    out_ids = os.path.join(work_dir, f"{prefix}.id_list")

    new_to_old_map = {}
    for line in id_mapping:
        old_id, new_id = line.strip().split('\t')
        new_to_old_map[new_id] = old_id

    with open(out_gff, 'w') as f_gff, \
         open(out_bed, 'w') as f_bed, \
         open(out_cds, 'w') as f_cds, \
         open(out_pep, 'w') as f_pep, \
         open(out_ids, 'w') as f_ids:

        for line in id_mapping:
            f_ids.write(line + "\n")

        f_gff.write("##gff-version 3\n")

        for gid, gene in tqdm(processed_genes.items(),
                              desc="         [Progress] Exporting   ",
                              unit=" gene", ncols=100, leave=False):
            scaf   = gene['scaf']
            strand = gene['strand']

            # Source field: original annotation source or '.' (default)
            gene_src = gene.get('source', '.') if keep_source else '.'

            g_attrs = f"ID={gid};Name={gid}"
            f_gff.write(
                f"{scaf}\t{gene_src}\tgene\t{gene['start']}\t{gene['end']}"
                f"\t.\t{strand}\t.\t{g_attrs}\n")

            for mid, mrna in gene['mrnas'].items():
                m_start = mrna['start']
                m_end   = mrna['end']

                m_attrs = f"ID={mid};Parent={gid};Name={mid}"
                # mRNA shares gene-level source (not stored separately)
                f_gff.write(
                    f"{scaf}\t{gene_src}\tmRNA\t{m_start}\t{m_end}"
                    f"\t.\t{strand}\t.\t{m_attrs}\n")

                # BED output
                if bed6:
                    # Standard BED6: chrom start end name score strand
                    f_bed.write(
                        f"{scaf}\t{m_start}\t{m_end}\t{mid}\t.\t{strand}\n")
                else:
                    # Extended format (original): includes old ID as col 7
                    old_mid = new_to_old_map.get(mid, ".")
                    f_bed.write(
                        f"{scaf}\t{m_start}\t{m_end}\t{mid}\t{strand}\t.\t{old_mid}\n")

                for sub in mrna['subfeats']:
                    s_type  = sub['type']
                    s_start = sub['start']
                    s_end   = sub['end']
                    s_phase = sub['phase'] if sub['phase'] is not None else '.'
                    s_src   = sub.get('source', '.') if keep_source else '.'

                    s_attrs = (f"ID={sub['attrs']['ID']};"
                               f"Parent={sub['attrs']['Parent']}")
                    f_gff.write(
                        f"{scaf}\t{s_src}\t{s_type}\t{s_start}\t{s_end}"
                        f"\t.\t{strand}\t{s_phase}\t{s_attrs}\n")

        for mid, seq in cds_records.items():
            f_cds.write(f">{mid}\n{seq}\n")

        for mid, seq in pep_records.items():
            f_pep.write(f">{mid}\n{seq}\n")

    if logger:
        logger.info(f"         ==> Generated {out_gff}")
        logger.info(f"         ==> Generated {out_bed}"
                    + (" (BED6 standard format)" if bed6
                       else " (extended format with original ID)"))
        logger.info(f"         ==> Generated {out_cds}")
        logger.info(f"         ==> Generated {out_pep}")
        logger.info(f"         ==> Generated {out_ids}")
