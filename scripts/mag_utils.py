"""Shared parsing utilities for MAG-module scripts."""


def parse_hmmsearch_tblout(tblout_path, evalue_threshold=1e-5):
    """Parse an hmmsearch --tblout file, returning one dict per passing hit.

    Each dict has keys: locus_tag, clb_gene, evalue, contig.
    contig is derived from the Prodigal protein name (contig_N → contig).
    Duplicate (locus_tag, clb_gene) pairs are deduplicated; the lowest
    e-value hit is kept.
    """
    best = {}  # (locus_tag, clb_gene) → lowest evalue hit dict
    with open(tblout_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 19:
                continue
            evalue = float(parts[4])
            if evalue > evalue_threshold:
                continue
            locus_tag = parts[0]
            clb_gene = parts[2]
            key = (locus_tag, clb_gene)
            if key not in best or evalue < best[key]["evalue"]:
                contig = "_".join(locus_tag.split("_")[:-1])
                best[key] = {
                    "locus_tag": locus_tag,
                    "clb_gene": clb_gene,
                    "evalue": evalue,
                    "contig": contig,
                }
    return list(best.values())
