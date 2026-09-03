#!/usr/bin/env python3
"""Build pks_mag_summary.tsv from per-sample MAG pipeline outputs."""

import argparse
import csv
import glob
import os
import sys
from pathlib import Path

from mag_utils import parse_hmmsearch_tblout

ENTEROBACTERALES = "o__Enterobacterales"

FIELDNAMES = [
    "sample", "bin_id", "taxonomy", "completeness", "contamination",
    "genome_size", "clb_genes_detected", "clb_genes", "best_evalue",
    "has_integrase", "has_transposase", "flanking_genes", "unexpected_taxon_flag",
]


def is_unexpected_taxon(taxonomy):
    """True only when taxonomy is classified but not Enterobacterales.
    Unclassified bins are unknown, not anomalous."""
    return taxonomy not in ("unclassified", "") and ENTEROBACTERALES not in taxonomy


def parse_checkm2(tsv_path):
    result = {}
    with open(tsv_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            result[row["Name"]] = {
                "completeness": float(row["Completeness"]),
                "contamination": float(row["Contamination"]),
                "genome_size": int(row.get("Genome_size", 0)),
            }
    return result


def parse_gtdbtk(summary_path):
    result = {}
    with open(summary_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            result[row["user_genome"]] = row.get("classification", "unclassified") or "unclassified"
    return result


def parse_context(context_path):
    has_int, has_tra, flank_set = False, False, set()
    with open(context_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            has_int |= row["has_integrase"].lower() == "true"
            has_tra |= row["has_transposase"].lower() == "true"
            if row.get("flanking_genes"):
                for g in row["flanking_genes"].split(";"):
                    if g.strip():
                        flank_set.add(g.strip())
    return has_int, has_tra, ";".join(sorted(flank_set))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkm2", required=True)
    parser.add_argument("--gtdbtk", required=True)
    parser.add_argument("--tblout_dir", required=True)
    parser.add_argument("--context_dir", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--evalue", type=float, default=1e-5)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    checkm2 = parse_checkm2(args.checkm2)
    gtdbtk = parse_gtdbtk(args.gtdbtk)

    rows = []
    for tblout_path in sorted(glob.glob(os.path.join(args.tblout_dir, "*.tblout"))):
        bin_id = Path(tblout_path).stem
        hits = parse_hmmsearch_tblout(tblout_path, args.evalue)
        if not hits:
            continue

        clb_genes_seen = {}
        for h in hits:
            gene = h["clb_gene"]
            clb_genes_seen[gene] = min(h["evalue"], clb_genes_seen.get(gene, float("inf")))
        clb_genes = sorted(clb_genes_seen.keys())
        best_evalue = min(clb_genes_seen.values())

        qc = checkm2.get(bin_id, {"completeness": 0.0, "contamination": 0.0, "genome_size": 0})
        taxonomy = gtdbtk.get(bin_id, "unclassified")
        unexpected = is_unexpected_taxon(taxonomy)

        context_path = os.path.join(args.context_dir, f"{bin_id}.context.tsv")
        if os.path.exists(context_path):
            has_int, has_tra, flanking = parse_context(context_path)
        else:
            has_int, has_tra, flanking = False, False, ""

        rows.append({
            "sample": args.sample,
            "bin_id": bin_id,
            "taxonomy": taxonomy,
            "completeness": qc["completeness"],
            "contamination": qc["contamination"],
            "genome_size": qc["genome_size"],
            "clb_genes_detected": len(clb_genes),
            "clb_genes": ",".join(clb_genes),
            "best_evalue": best_evalue,
            "has_integrase": has_int,
            "has_transposase": has_tra,
            "flanking_genes": flanking,
            "unexpected_taxon_flag": unexpected,
        })

    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDNAMES, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    if not rows:
        print(f"No pks+ bins found for {args.sample}", file=sys.stderr)
    else:
        print(f"Wrote {len(rows)} pks+ bins to {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
