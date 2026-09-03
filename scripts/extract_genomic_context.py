#!/usr/bin/env python3
"""Extract genomic context around clb gene hits in a pks+ MAG bin."""

import argparse
import csv
import re
import sys


def parse_prokka_gff(gff_path):
    genes = []
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] not in ("gene", "CDS"):
                continue
            attr = dict(re.findall(r'(\w+)=([^;]+)', parts[8]))
            genes.append({
                "contig": parts[0],
                "start": int(parts[3]),
                "end": int(parts[4]),
                "strand": parts[6],
                "locus_tag": attr.get("ID", attr.get("locus_tag", "")),
                "gene": attr.get("gene", ""),
                "product": attr.get("product", ""),
            })
    return genes


def parse_hmmsearch_tblout(tblout_path, evalue_threshold=1e-5):
    hits = []
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
            target = parts[0]
            query = parts[2]
            contig = "_".join(target.split("_")[:-1])
            hits.append({
                "locus_tag": target,
                "clb_gene": query,
                "evalue": evalue,
                "contig": contig,
            })
    return hits


def get_flanking(genes, contig, pos, window=50000):
    return [
        g for g in genes
        if g["contig"] == contig
        and not (g["end"] < pos - window or g["start"] > pos + window)
    ]


def is_integrase(gene):
    keywords = {"integrase", "recombinase", "resolvase"}
    text = (gene["gene"] + " " + gene["product"]).lower()
    return any(kw in text for kw in keywords)


def is_transposase(gene):
    keywords = {"transposase", "insertion element", "is element"}
    text = (gene["gene"] + " " + gene["product"]).lower()
    return any(kw in text for kw in keywords)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gff", required=True)
    parser.add_argument("--tblout", required=True)
    parser.add_argument("--evalue", type=float, default=1e-5)
    parser.add_argument("--window", type=int, default=50000)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    genes = parse_prokka_gff(args.gff)
    hits = parse_hmmsearch_tblout(args.tblout, args.evalue)

    tag_to_gene = {g["locus_tag"]: g for g in genes}

    rows = []
    for hit in hits:
        anchor = tag_to_gene.get(hit["locus_tag"])
        mid = (anchor["start"] + anchor["end"]) // 2 if anchor else 0
        flanking = get_flanking(genes, hit["contig"], mid, args.window)
        has_int = any(is_integrase(g) for g in flanking)
        has_tra = any(is_transposase(g) for g in flanking)
        flank_names = ";".join(
            g["gene"] or g["product"]
            for g in flanking
            if g["locus_tag"] != hit["locus_tag"]
        )
        rows.append({
            "locus_tag": hit["locus_tag"],
            "clb_gene": hit["clb_gene"],
            "evalue": hit["evalue"],
            "contig": hit["contig"],
            "has_integrase": has_int,
            "has_transposase": has_tra,
            "flanking_genes": flank_names,
        })

    fieldnames = ["locus_tag", "clb_gene", "evalue", "contig",
                  "has_integrase", "has_transposase", "flanking_genes"]
    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} context rows to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
