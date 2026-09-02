# MAG Module Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement `Modules/pks_mag.nf` and the MAG summary table process — an optional metagenome assembly → binning → pks protein HMM detection branch for pksProfiler.

**Architecture:** nf-core DSL2 modules imported into a `pksMAG` wrapper workflow in `Modules/pks_mag.nf`. Custom processes handle genomic context extraction and MAG summary table generation. The branch is triggered by `params.enable_mags = true` in `main.nf` and runs on `MAPPED_READS` (all clean non-host reads, not entero-filtered).

**Tech Stack:** Nextflow DSL2, nf-core modules (MEGAHIT, minimap2, samtools, MetaBAT2, CheckM2, GTDB-Tk, Prodigal, HMMER, Prokka), Python 3.10, MAFFT, HMMER3.

---

## File Structure

| File | Action | Responsibility |
|------|--------|----------------|
| `ref/hmm/clb_all_protein.hmm` | **Build** | colibactin protein HMM profiles (19 genes, multi-strain) |
| `modules/nf-core/*/main.nf` | **Install** | versioned nf-core tool wrappers |
| `conf/base.config` | **Modify** | 4 new MAG resource labels |
| `scripts/extract_genomic_context.py` | **Create** | parse Prokka GFF + hmmsearch tblout → per-bin context TSV |
| `scripts/build_mag_summary.py` | **Create** | aggregate CheckM2 + GTDB-Tk + hmmsearch + context → `pks_mag_summary.tsv` |
| `Modules/pks_mag.nf` | **Create** | `pksMAG` workflow + `extractGenomicContext` + `magSummaryTable` processes |
| `main.nf` | **Modify** | `params.enable_mags`, `include`, MAG branch wiring |

---

## Task 1: Build colibactin protein HMM (offline prerequisite)

**Files:**
- Create: `ref/hmm/clb_all_protein.hmm` + `.h3{f,i,m,p}`
- Create: `scripts/build_clb_protein_hmm.sh` (reproducibility script)

This is an offline one-time step. The resulting HMM file is committed to the repo so pipeline users never need to rebuild it.

- [ ] **Step 1: Set up build environment**

```bash
mamba create -n hmm_build -c bioconda -c conda-forge \
    hmmer=3.3.2 mafft=7.526 entrez-direct=22.4
conda activate hmm_build
```

- [ ] **Step 2: Fetch clbA–clbS proteins from NCBI for multiple strains**

```bash
mkdir -p /tmp/clb_seqs && cd /tmp/clb_seqs

# Strain accession list (pks+ E. coli and Klebsiella)
# IHE3034 (NC_017628.1), CFT073 (AE014075.1), Nissle 1917 (CP007799.1),
# 536 (CP000247.1), APEC O1 (CP000468.1)
STRAINS=(
  "IHE3034[Organism] AND pks island[Title]"
  "CFT073[Organism]"
  "Nissle 1917[Organism]"
)

for gene in A B C D E F G H I J K L M N O P Q R S; do
  echo "Fetching clb${gene}..."
  esearch -db protein \
    -query "clb${gene}[Gene Name] Escherichia coli[Organism]" |
    efetch -format fasta >> "clb${gene}.fasta"

  # Also fetch Klebsiella pneumoniae homologs
  esearch -db protein \
    -query "clb${gene}[Gene Name] Klebsiella[Organism]" |
    efetch -format fasta >> "clb${gene}.fasta"

  # Remove empty sequences and deduplicate
  seqkit rmdup -s "clb${gene}.fasta" -o "clb${gene}.dedup.fasta"
done
```

Expected: each `clbX.dedup.fasta` has 3–10 sequences.

- [ ] **Step 3: Align and build per-gene HMMs**

```bash
cd /tmp/clb_seqs
mkdir -p hmms

for gene in A B C D E F G H I J K L M N O P Q R S; do
  if [[ ! -s "clb${gene}.dedup.fasta" ]]; then
    echo "WARNING: no sequences for clb${gene}, skipping"
    continue
  fi

  # Multiple sequence alignment
  mafft --auto --thread 4 "clb${gene}.dedup.fasta" > "clb${gene}.msa"

  # Build HMM profile
  hmmbuild \
    --cpu 4 \
    -n "clb${gene}" \
    "hmms/clb${gene}.hmm" \
    "clb${gene}.msa"
done
```

- [ ] **Step 4: Concatenate and press**

```bash
cd /tmp/clb_seqs
cat hmms/clb*.hmm > clb_all_protein.hmm
hmmpress clb_all_protein.hmm
ls -lh clb_all_protein.hmm clb_all_protein.hmm.h3{f,i,m,p}
```

Expected output (5 files): `clb_all_protein.hmm`, `.h3f`, `.h3i`, `.h3m`, `.h3p`

- [ ] **Step 5: Validate with a test search**

```bash
# Test against IHE3034 proteins (should hit all 19 clb genes)
# Download IHE3034 proteins: NC_017628.1 → fetch protein FASTA
esearch -db protein -query "NC_017628.1[Accession]" | efetch -format fasta > ihe3034.faa

hmmsearch \
  --tblout test_hits.tbl \
  --noali \
  -E 1e-5 \
  clb_all_protein.hmm \
  ihe3034.faa

grep -v '^#' test_hits.tbl | awk '{print $3}' | sort -u
```

Expected: all 19 clbA–clbS names appear in output.

- [ ] **Step 6: Copy to repo and commit**

```bash
cd /tscc/lustre/restricted/alexandrov-ddn/users/l1joseph/pksProfiler
cp /tmp/clb_seqs/clb_all_protein.hmm ref/hmm/
cp /tmp/clb_seqs/clb_all_protein.hmm.h3{f,i,m,p} ref/hmm/

# Save build script for reproducibility
cp /tmp/clb_seqs build_clb_protein_hmm.sh  # (copy commands from steps 2-4)

git add ref/hmm/clb_all_protein.hmm ref/hmm/clb_all_protein.hmm.h3{f,i,m,p}
git commit -m "Add colibactin protein HMM (19 genes, multi-strain)"
```

---

## Task 2: Add MAG resource labels to conf/base.config

**Files:**
- Modify: `conf/base.config`

- [ ] **Step 1: Open conf/base.config and locate the closing `}`**

The file ends with the `entero_filter` label block. Add four new labels before the final `}`.

- [ ] **Step 2: Add resource labels**

In `conf/base.config`, immediately before the final `}`:

```groovy
    withLabel:mag_assembly {
        cpus   = 16
        memory = { 64.GB  * task.attempt }
        time   = { 12.h   * task.attempt }
    }
    withLabel:mag_binning {
        cpus   = 8
        memory = { 32.GB  * task.attempt }
        time   = { 4.h    * task.attempt }
    }
    withLabel:mag_gtdbtk {
        cpus   = 16
        memory = { 200.GB * task.attempt }
        time   = { 24.h   * task.attempt }
    }
    withLabel:mag_hmm {
        cpus   = 8
        memory = { 16.GB  * task.attempt }
        time   = { 4.h    * task.attempt }
    }
```

- [ ] **Step 3: Verify syntax**

```bash
nextflow config -show-profiles 2>&1 | head -5
```

Expected: no Groovy parse errors.

- [ ] **Step 4: Commit**

```bash
git add conf/base.config
git commit -m "Add MAG resource labels to base.config"
```

---

## Task 3: Install nf-core modules

**Files:**
- Create: `modules/nf-core/*/main.nf` (11 modules, via CLI)

nf-core modules install into `modules/nf-core/<tool>/main.nf`. Each module carries its own versioned conda spec — no new `conda_envs/` YML needed.

- [ ] **Step 1: Install nf-core CLI if not present**

```bash
pip install nf-core
nf-core --version
```

Expected: `nf-core, version 3.x.x`

- [ ] **Step 2: Install all 11 modules**

Run from the project root (`/tscc/lustre/restricted/alexandrov-ddn/users/l1joseph/pksProfiler`):

```bash
nf-core modules install megahit
nf-core modules install minimap2/align
nf-core modules install samtools/sort
nf-core modules install samtools/index
nf-core modules install metabat2/jgisummarizebamcontigdepths
nf-core modules install metabat2/metabat2
nf-core modules install checkm2/predict
nf-core modules install gtdbtk/classifywf
nf-core modules install prodigal
nf-core modules install hmmer/hmmsearch
nf-core modules install prokka
```

- [ ] **Step 3: Verify installation**

```bash
ls modules/nf-core/
# Expected directories: megahit  minimap2  samtools  metabat2  checkm2  gtdbtk  prodigal  hmmer  prokka

# Check each module has main.nf
find modules/nf-core -name "main.nf" | sort
```

- [ ] **Step 4: Read each module's input/output signature**

This is critical — nf-core module interfaces vary by version. Read and note the exact input tuple structure for each module before writing `pks_mag.nf`:

```bash
grep -A 10 "^    input:" modules/nf-core/megahit/main.nf
grep -A 10 "^    input:" modules/nf-core/minimap2/align/main.nf
grep -A 10 "^    input:" modules/nf-core/metabat2/metabat2/main.nf
grep -A 10 "^    input:" modules/nf-core/checkm2/predict/main.nf
grep -A 10 "^    input:" modules/nf-core/gtdbtk/classifywf/main.nf
grep -A 10 "^    input:" modules/nf-core/prodigal/main.nf
grep -A 10 "^    input:" modules/nf-core/hmmer/hmmsearch/main.nf
grep -A 10 "^    input:" modules/nf-core/prokka/main.nf
```

- [ ] **Step 5: Commit installed modules**

```bash
git add modules/
git commit -m "Install nf-core modules for MAG branch"
```

---

## Task 4: Create scripts/extract_genomic_context.py

**Files:**
- Create: `scripts/extract_genomic_context.py`

Parses Prokka GFF annotation and hmmsearch tblout for a single bin. For each clb protein hit, extracts gene names in the ±50 kb flanking region and checks for integrase/transposase annotations (using Prokka's gene name annotations, which annotate these from its internal databases).

- [ ] **Step 1: Write the script**

```python
#!/usr/bin/env python3
"""Extract genomic context around clb hits in a MAG bin.

For each clb protein hit in the hmmsearch tblout, reports genes in the
±50 kb flanking region from the Prokka GFF annotation, and flags whether
the region contains predicted integrases or transposases.
"""

import argparse
import sys
from collections import defaultdict
from pathlib import Path


INTEGRASE_KEYWORDS = {"integrase", "int", "inti", "intI", "site-specific recombinase"}
TRANSPOSASE_KEYWORDS = {"transposase", "tnp", "tnpA", "tnpB", "IS element"}
CLB_EVALUE = 1e-5


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--gff", required=True, type=Path,
                   help="Prokka GFF3 annotation for this bin")
    p.add_argument("--tblout", required=True, type=Path,
                   help="hmmsearch --tblout output against clb_all_protein.hmm")
    p.add_argument("--flank-bp", type=int, default=50_000)
    p.add_argument("--output", required=True, type=Path)
    return p.parse_args()


def parse_prokka_gff(gff_path):
    """Return list of (contig, start, end, gene_name, product) for CDS features."""
    features = []
    with gff_path.open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            contig = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            attrs = {}
            for a in parts[8].split(";"):
                if "=" in a:
                    k, v = a.split("=", 1)
                    attrs[k.strip()] = v.strip()
            gene_name = attrs.get("gene", "")
            product = attrs.get("product", "")
            locus = attrs.get("ID", "")
            features.append((contig, start, end, gene_name, product, locus))
    return features


def parse_hmmsearch_tblout(tblout_path, evalue_thresh=CLB_EVALUE):
    """Return dict: locus_tag → (clb_gene, evalue) for hits passing threshold."""
    hits = {}
    with tblout_path.open() as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if len(parts) < 10:
                continue
            target = parts[0]   # protein locus tag
            hmm_name = parts[2] # e.g. clbA, clbB
            evalue = float(parts[4])
            if evalue <= evalue_thresh:
                if target not in hits or evalue < hits[target][1]:
                    hits[target] = (hmm_name, evalue)
    return hits


def get_contig_for_locus(features, locus_tag):
    """Find which contig a locus tag belongs to."""
    for contig, start, end, gene, product, locus in features:
        if locus == locus_tag:
            return contig, start, end
    return None, None, None


def get_flanking(features, contig, hit_start, hit_end, flank_bp):
    """Return list of (gene_name, product) in flank_bp around hit region."""
    region_start = max(1, hit_start - flank_bp)
    region_end = hit_end + flank_bp
    return [
        (gene, product)
        for c, s, e, gene, product, locus in features
        if c == contig and s >= region_start and e <= region_end
    ]


def is_integrase(gene_name, product):
    name = (gene_name + " " + product).lower()
    return any(kw.lower() in name for kw in INTEGRASE_KEYWORDS)


def is_transposase(gene_name, product):
    name = (gene_name + " " + product).lower()
    return any(kw.lower() in name for kw in TRANSPOSASE_KEYWORDS)


def main():
    args = parse_args()

    features = parse_prokka_gff(args.gff)
    hits = parse_hmmsearch_tblout(args.tblout)

    rows = []
    for locus_tag, (clb_gene, evalue) in hits.items():
        contig, hit_start, hit_end = get_contig_for_locus(features, locus_tag)
        if contig is None:
            contig, hit_start, hit_end = "unknown", 0, 0

        flanking = get_flanking(features, contig, hit_start, hit_end, args.flank_bp)

        has_integrase = any(is_integrase(g, p) for g, p in flanking)
        has_transposase = any(is_transposase(g, p) for g, p in flanking)

        flanking_gene_names = [g for g, p in flanking if g][:30]

        rows.append({
            "locus_tag": locus_tag,
            "clb_gene": clb_gene,
            "evalue": evalue,
            "contig": contig,
            "has_integrase": has_integrase,
            "has_transposase": has_transposase,
            "flanking_genes": ";".join(flanking_gene_names),
        })

    with args.output.open("w") as fh:
        fh.write("locus_tag\tclb_gene\tevalue\tcontig\t"
                 "has_integrase\thas_transposase\tflanking_genes\n")
        for row in rows:
            fh.write(
                f"{row['locus_tag']}\t{row['clb_gene']}\t{row['evalue']}\t"
                f"{row['contig']}\t{row['has_integrase']}\t"
                f"{row['has_transposase']}\t{row['flanking_genes']}\n"
            )

    if not rows:
        sys.exit(0)


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Write unit tests**

Create `tests/test_extract_genomic_context.py`:

```python
import textwrap
from pathlib import Path
import pytest
import sys
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))
from extract_genomic_context import (
    parse_prokka_gff, parse_hmmsearch_tblout,
    get_flanking, is_integrase, is_transposase,
)


GFF_CONTENT = textwrap.dedent("""\
    ##gff-version 3
    bin.1\tProdigal:002006483\tCDS\t100\t900\t.\t+\t0\tID=PROKKA_00001;gene=clbA;product=ClbA
    bin.1\tProdigal:002006483\tCDS\t1000\t1800\t.\t+\t0\tID=PROKKA_00002;gene=intI;product=integrase
    bin.1\tProdigal:002006483\tCDS\t200000\t201000\t.\t-\t0\tID=PROKKA_00003;gene=tnpA;product=transposase
""")

TBLOUT_CONTENT = textwrap.dedent("""\
    # hmmsearch
    PROKKA_00001  -  clbA  -  1.2e-10  50.0  0.1  1  0  0  1  1  1  1  -
""")


def test_parse_gff(tmp_path):
    gff = tmp_path / "test.gff"
    gff.write_text(GFF_CONTENT)
    features = parse_prokka_gff(gff)
    assert len(features) == 3
    assert features[0][0] == "bin.1"   # contig
    assert features[0][3] == "clbA"    # gene name


def test_parse_tblout(tmp_path):
    tbl = tmp_path / "test.tbl"
    tbl.write_text(TBLOUT_CONTENT)
    hits = parse_hmmsearch_tblout(tbl, evalue_thresh=1e-5)
    assert "PROKKA_00001" in hits
    assert hits["PROKKA_00001"][0] == "clbA"


def test_flanking_integrase(tmp_path):
    gff = tmp_path / "test.gff"
    gff.write_text(GFF_CONTENT)
    features = parse_prokka_gff(gff)
    # clbA hit at 100-900, integrase at 1000-1800 — within 50kb
    flanking = get_flanking(features, "bin.1", 100, 900, 50_000)
    assert any(is_integrase(g, p) for g, p in flanking)


def test_flanking_transposase_out_of_range(tmp_path):
    gff = tmp_path / "test.gff"
    gff.write_text(GFF_CONTENT)
    features = parse_prokka_gff(gff)
    # transposase at 200000 is >50kb from clbA hit at 100-900
    flanking = get_flanking(features, "bin.1", 100, 900, 50_000)
    assert not any(is_transposase(g, p) for g, p in flanking)
```

- [ ] **Step 3: Run tests**

```bash
cd /tscc/lustre/restricted/alexandrov-ddn/users/l1joseph/pksProfiler
conda activate blahhrd  # any env with pytest + python3
pytest tests/test_extract_genomic_context.py -v
```

Expected: 4 tests pass.

- [ ] **Step 4: Commit**

```bash
git add scripts/extract_genomic_context.py tests/test_extract_genomic_context.py
git commit -m "Add extract_genomic_context.py for MAG branch"
```

---

## Task 5: Create scripts/build_mag_summary.py

**Files:**
- Create: `scripts/build_mag_summary.py`

Aggregates per-sample CheckM2 reports, GTDB-Tk classification, hmmsearch tblouts, and context TSVs into one cohort-level `pks_mag_summary.tsv`. Only bins with ≥1 clb gene hit appear in the output.

- [ ] **Step 1: Write the script**

```python
#!/usr/bin/env python3
"""Build cohort-level pks_mag_summary.tsv from per-sample MAG outputs.

Input files (all passed as lists via CLI):
  --checkm2   : CheckM2 quality_report.tsv files (one per sample)
  --gtdbtk    : GTDB-Tk gtdbtk.bac120.summary.tsv files (one per sample)
  --hmmsearch : hmmsearch tblout files (one per bin)
  --context   : extract_genomic_context output TSVs (one per pks+ bin)
  --output    : path for pks_mag_summary.tsv
"""

import argparse
import csv
from collections import defaultdict
from pathlib import Path


CLB_EVALUE_THRESH = 1e-5
ENTEROBACTERALES_KEYWORDS = {
    "Enterobacteriaceae", "Enterobacterales", "Escherichia",
    "Klebsiella", "Salmonella", "Shigella", "Yersinia",
    "Serratia", "Enterobacter", "Citrobacter", "Cronobacter",
}


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--checkm2",   nargs="+", required=True, type=Path)
    p.add_argument("--gtdbtk",    nargs="+", required=True, type=Path)
    p.add_argument("--hmmsearch", nargs="+", required=True, type=Path)
    p.add_argument("--context",   nargs="+", required=False, default=[], type=Path)
    p.add_argument("--output",    required=True, type=Path)
    return p.parse_args()


def parse_checkm2(paths):
    """Return dict: bin_id → {completeness, contamination, genome_size}."""
    result = {}
    for path in paths:
        with path.open(newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                bin_id = row["Name"]
                result[bin_id] = {
                    "completeness": float(row.get("Completeness", 0)),
                    "contamination": float(row.get("Contamination", 0)),
                    "genome_size": int(row.get("Genome_Size", 0)),
                }
    return result


def parse_gtdbtk(paths):
    """Return dict: bin_id → taxonomy string."""
    result = {}
    for path in paths:
        if not path.exists() or path.stat().st_size == 0:
            continue
        with path.open(newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                bin_id = row.get("user_genome", row.get("Name", ""))
                taxonomy = row.get("classification", "")
                result[bin_id] = taxonomy
    return result


def parse_hmmsearch_tblout(paths, evalue_thresh=CLB_EVALUE_THRESH):
    """Return dict: bin_id → {clb_genes: set, best_evalue: float}.

    Tblout file names are expected to be <sample>.<bin_id>.tblout
    so bin_id is extracted from the filename stem.
    """
    result = defaultdict(lambda: {"clb_genes": set(), "best_evalue": float("inf")})
    for path in paths:
        # Filename: sampleID.bin_id.tblout  →  bin_id is stem minus last extension
        bin_id = path.stem  # e.g. CM017.bin.1
        with path.open() as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split()
                if len(parts) < 10:
                    continue
                hmm_name = parts[2]
                evalue = float(parts[4])
                if evalue <= evalue_thresh:
                    result[bin_id]["clb_genes"].add(hmm_name)
                    if evalue < result[bin_id]["best_evalue"]:
                        result[bin_id]["best_evalue"] = evalue
    return dict(result)


def parse_context(paths):
    """Return dict: bin_id → {has_integrase, has_transposase, flanking_genes}."""
    result = {}
    for path in paths:
        bin_id = path.stem  # sampleID.bin_id.context
        has_integrase = False
        has_transposase = False
        all_flanking = set()
        with path.open(newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                if row.get("has_integrase", "False") == "True":
                    has_integrase = True
                if row.get("has_transposase", "False") == "True":
                    has_transposase = True
                for gene in row.get("flanking_genes", "").split(";"):
                    if gene.strip():
                        all_flanking.add(gene.strip())
        result[bin_id] = {
            "has_integrase": has_integrase,
            "has_transposase": has_transposase,
            "flanking_genes": ";".join(sorted(all_flanking)[:30]),
        }
    return result


def is_unexpected_taxon(taxonomy):
    """True if GTDB-Tk assigns outside Enterobacterales."""
    if not taxonomy:
        return False
    return not any(kw in taxonomy for kw in ENTEROBACTERALES_KEYWORDS)


COLUMNS = [
    "sample", "bin_id", "taxonomy", "completeness", "contamination",
    "genome_size", "clb_genes_detected", "clb_genes", "best_evalue",
    "has_integrase", "has_transposase", "flanking_genes", "unexpected_taxon_flag",
]


def main():
    args = parse_args()

    checkm2 = parse_checkm2(args.checkm2)
    gtdbtk = parse_gtdbtk(args.gtdbtk)
    hmm = parse_hmmsearch_tblout(args.hmmsearch)
    context = parse_context(args.context) if args.context else {}

    # Only emit bins with ≥1 clb hit
    pks_bins = {bin_id: data for bin_id, data in hmm.items()
                if data["clb_genes"]}

    rows = []
    for bin_id, hdata in pks_bins.items():
        # bin_id format: sampleID.bin_X → split on first dot for sample
        parts = bin_id.split(".", 1)
        sample = parts[0] if len(parts) == 2 else bin_id
        inner_bin = parts[1] if len(parts) == 2 else bin_id

        qc = checkm2.get(inner_bin, {"completeness": "", "contamination": "", "genome_size": ""})
        taxonomy = gtdbtk.get(inner_bin, "")
        ctx = context.get(bin_id, {"has_integrase": "", "has_transposase": "", "flanking_genes": ""})

        clb_genes_sorted = ",".join(sorted(hdata["clb_genes"]))
        best_ev = "" if hdata["best_evalue"] == float("inf") else f"{hdata['best_evalue']:.2e}"

        rows.append({
            "sample": sample,
            "bin_id": inner_bin,
            "taxonomy": taxonomy,
            "completeness": qc["completeness"],
            "contamination": qc["contamination"],
            "genome_size": qc["genome_size"],
            "clb_genes_detected": len(hdata["clb_genes"]),
            "clb_genes": clb_genes_sorted,
            "best_evalue": best_ev,
            "has_integrase": ctx["has_integrase"],
            "has_transposase": ctx["has_transposase"],
            "flanking_genes": ctx["flanking_genes"],
            "unexpected_taxon_flag": is_unexpected_taxon(taxonomy),
        })

    rows.sort(key=lambda r: (r["sample"], -int(r["clb_genes_detected"])))

    with args.output.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} pks+ MAG rows to {args.output}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Write tests**

Create `tests/test_build_mag_summary.py`:

```python
import textwrap
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))
from build_mag_summary import (
    parse_checkm2, parse_gtdbtk, parse_hmmsearch_tblout,
    parse_context, is_unexpected_taxon,
)


CHECKM2_TSV = textwrap.dedent("""\
    Name\tCompleteness\tcontamination\tGenome_Size
    bin.1\t95.2\t1.3\t5200000
    bin.2\t42.0\t5.1\t1800000
""")

GTDBTK_TSV = textwrap.dedent("""\
    user_genome\tclassification
    bin.1\td__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli
    bin.2\td__Bacteria;p__Firmicutes;c__Bacilli;o__Lactobacillales;f__Lactobacillaceae;g__Lactobacillus;s__Lactobacillus acidophilus
""")

TBLOUT = textwrap.dedent("""\
    # hmmsearch
    PROKKA_00001  -  clbA  -  1.2e-10  50.0  0.1  1  0  0  1  1  1  1  -
    PROKKA_00002  -  clbB  -  3.4e-15  80.2  0.0  1  0  0  1  1  1  1  -
""")

CONTEXT_TSV = textwrap.dedent("""\
    locus_tag\tclb_gene\tevalue\tcontig\thas_integrase\thas_transposase\tflanking_genes
    PROKKA_00001\tclbA\t1.2e-10\tbin.1_contig_1\tTrue\tFalse\tint1;recA
""")


def test_parse_checkm2(tmp_path):
    f = tmp_path / "quality_report.tsv"
    f.write_text(CHECKM2_TSV)
    data = parse_checkm2([f])
    assert data["bin.1"]["completeness"] == 95.2
    assert "bin.2" in data


def test_parse_gtdbtk(tmp_path):
    f = tmp_path / "gtdbtk.bac120.summary.tsv"
    f.write_text(GTDBTK_TSV)
    data = parse_gtdbtk([f])
    assert "Escherichia" in data["bin.1"]
    assert "Lactobacillus" in data["bin.2"]


def test_parse_hmmsearch(tmp_path):
    f = tmp_path / "CM017.bin.1.tblout"
    f.write_text(TBLOUT)
    data = parse_hmmsearch_tblout([f])
    assert "clbA" in data["CM017.bin.1"]["clb_genes"]
    assert "clbB" in data["CM017.bin.1"]["clb_genes"]
    assert data["CM017.bin.1"]["best_evalue"] < 1e-5


def test_unexpected_taxon():
    assert not is_unexpected_taxon("o__Enterobacterales;f__Enterobacteriaceae")
    assert is_unexpected_taxon("o__Lactobacillales;f__Lactobacillaceae")
    assert not is_unexpected_taxon("")
```

- [ ] **Step 3: Run tests**

```bash
pytest tests/test_build_mag_summary.py -v
```

Expected: 4 tests pass.

- [ ] **Step 4: Commit**

```bash
git add scripts/build_mag_summary.py tests/test_build_mag_summary.py
git commit -m "Add build_mag_summary.py for MAG summary table"
```

---

## Task 6: Create Modules/pks_mag.nf — custom processes

**Files:**
- Create: `Modules/pks_mag.nf` (custom processes only in this task)

Write the two custom processes (`extractGenomicContext` and `magSummaryTable`) that wrap the scripts from Tasks 4 and 5. Follow pksProfiler's existing process style: 2-space indent in script block, `scratch true`, `publishDir`, `set -euo pipefail`.

- [ ] **Step 1: Create Modules/pks_mag.nf with custom processes**

```groovy
nextflow.enable.dsl = 2

// ---------- nf-core module includes ----------
// (added in Task 7)

// ---------- Custom processes ----------

process extractGenomicContext {
  label 'mag_hmm'
  scratch true
  publishDir "${params.outdir}/pks_summary/mag/${sampleID}", mode: 'copy'
  conda 'conda-forge::python=3.10'

  input:
  tuple val(sampleID), val(binID), path(gff), path(tblout)

  output:
  tuple val(sampleID), val(binID), path("${sampleID}.${binID}.context.tsv")

  script:
  """
  set -euo pipefail

  python "${params.scripts}/extract_genomic_context.py" \
    --gff "${gff}" \
    --tblout "${tblout}" \
    --flank-bp 50000 \
    --output "${sampleID}.${binID}.context.tsv"
  """
}


process magSummaryTable {
  label 'mag_hmm'
  scratch true
  publishDir "${params.outdir}/pks_summary/mag", mode: 'copy'
  conda 'conda-forge::python=3.10'

  input:
  path checkm2_reports
  path gtdbtk_reports
  path hmmsearch_tblouts
  path context_tsvs

  output:
  path "pks_mag_summary.tsv"

  script:
  def checkm2_inputs  = checkm2_reports instanceof List  ? checkm2_reports.join(' ')  : checkm2_reports
  def gtdbtk_inputs   = gtdbtk_reports instanceof List   ? gtdbtk_reports.join(' ')   : gtdbtk_reports
  def hmm_inputs      = hmmsearch_tblouts instanceof List ? hmmsearch_tblouts.join(' ') : hmmsearch_tblouts
  def context_inputs  = context_tsvs instanceof List     ? context_tsvs.join(' ')     : context_tsvs

  """
  set -euo pipefail

  python "${params.scripts}/build_mag_summary.py" \
    --checkm2  ${checkm2_inputs} \
    --gtdbtk   ${gtdbtk_inputs} \
    --hmmsearch ${hmm_inputs} \
    --context  ${context_inputs} \
    --output pks_mag_summary.tsv
  """
}
```

- [ ] **Step 2: Verify the file parses**

```bash
nextflow run /dev/null --preview 2>&1 | head -3
# Or just check for syntax errors:
groovy -e 'new File("Modules/pks_mag.nf").text' 2>&1 | head -5
```

Expected: no Groovy parse errors.

- [ ] **Step 3: Commit (partial file — workflow added next task)**

```bash
git add Modules/pks_mag.nf
git commit -m "Add custom processes to pks_mag.nf"
```

---

## Task 7: Complete Modules/pks_mag.nf — pksMAG workflow

**Files:**
- Modify: `Modules/pks_mag.nf` (add nf-core includes + pksMAG workflow)

**Before starting:** Read the installed nf-core module signatures from Task 3, Step 4. The input tuple structures below match the standard nf-core convention — verify they match your installed versions and adjust if they differ.

nf-core modules use `tuple val(meta), path(...)` where `meta` is a Groovy map (`[id: sampleID]`). Our pipeline uses `tuple val(sampleID), path(...)` — a one-liner adapter converts between them at the workflow boundary.

- [ ] **Step 1: Add nf-core includes at the top of Modules/pks_mag.nf**

Insert after the `nextflow.enable.dsl = 2` line, before the custom processes:

```groovy
// ---------- nf-core module includes ----------
include { MEGAHIT                              } from '../modules/nf-core/megahit/main'
include { MINIMAP2_ALIGN                       } from '../modules/nf-core/minimap2/align/main'
include { SAMTOOLS_SORT                        } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX                       } from '../modules/nf-core/samtools/index/main'
include { METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS } from '../modules/nf-core/metabat2/jgisummarizebamcontigdepths/main'
include { METABAT2_METABAT2                    } from '../modules/nf-core/metabat2/metabat2/main'
include { CHECKM2_PREDICT                      } from '../modules/nf-core/checkm2/predict/main'
include { GTDBTK_CLASSIFYWF                    } from '../modules/nf-core/gtdbtk/classifywf/main'
include { PRODIGAL                             } from '../modules/nf-core/prodigal/main'
include { HMMER_HMMSEARCH                      } from '../modules/nf-core/hmmer/hmmsearch/main'
include { PROKKA                               } from '../modules/nf-core/prokka/main'
```

- [ ] **Step 2: Add the pksMAG workflow at the end of Modules/pks_mag.nf**

```groovy
workflow pksMAG {
  take:
  reads  // tuple val(sampleID), path(fastq_gz)

  main:
  // ── Channel adapter: sampleID string → nf-core meta map ──────────────────
  reads
    .map { sampleID, fastq ->
      tuple([id: sampleID, single_end: true], [fastq])
    }
    .set { reads_meta }

  // ── Assembly ─────────────────────────────────────────────────────────────
  MEGAHIT(reads_meta)

  // ── Map reads back to contigs (for binning depth) ─────────────────────
  // MINIMAP2_ALIGN: (meta, reads), reference_path, bam_format, cigar_paf_format, cigar_bam
  reads_meta
    .join(MEGAHIT.out.contigs)
    .map { meta, reads, contigs -> tuple(meta, reads, contigs, true, false, false) }
    .set { minimap2_input }

  MINIMAP2_ALIGN(minimap2_input.map { m, r, ref, b, c1, c2 -> tuple(m, r) },
                 minimap2_input.map { m, r, ref, b, c1, c2 -> ref },
                 true, false, false)

  SAMTOOLS_SORT(MINIMAP2_ALIGN.out.bam, [[], []])
  SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)

  // ── Contig depth for binning ───────────────────────────────────────────
  SAMTOOLS_SORT.out.bam
    .join(SAMTOOLS_INDEX.out.bai)
    .set { bam_bai }

  METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS(bam_bai)

  // ── Binning ───────────────────────────────────────────────────────────
  MEGAHIT.out.contigs
    .join(METABAT2_JGISUMMARIZEBAMCONTIGDEPTHS.out.depth)
    .set { binning_input }

  METABAT2_METABAT2(binning_input)

  // ── QC and taxonomy (operate on all bins per sample) ──────────────────
  def checkm2_db = file(params.checkm2_db, checkIfExists: true)
  def gtdbtk_db  = file(params.gtdbtk_db,  checkIfExists: true)

  CHECKM2_PREDICT(METABAT2_METABAT2.out.fasta, checkm2_db)
  GTDBTK_CLASSIFYWF(METABAT2_METABAT2.out.fasta, gtdbtk_db)

  // ── Protein prediction on each bin (flatMap bins → one channel per bin) ──
  def clb_protein_hmm = file(params.clb_protein_hmm, checkIfExists: true)

  METABAT2_METABAT2.out.fasta
    .flatMap { meta, bins ->
      (bins instanceof List ? bins : [bins]).collect { bin ->
        tuple(meta, bin)
      }
    }
    .set { per_bin }

  PRODIGAL(per_bin, 'gff')

  // ── clb protein HMM search (all bins, no taxonomy pre-filter) ─────────
  PRODIGAL.out.amino_acid_fasta
    .map { meta, faa -> tuple(meta, clb_protein_hmm, faa, []) }
    .set { hmmsearch_input }

  HMMER_HMMSEARCH(hmmsearch_input)

  // ── Filter to pks+ bins (≥1 clb hit) ─────────────────────────────────
  HMMER_HMMSEARCH.out.target_summary
    .filter { meta, tblout ->
      tblout.text.readLines().any { line ->
        !line.startsWith('#') && line.split(/\s+/).size() >= 5 &&
        line.split(/\s+/)[4].toDouble() <= 1e-5
      }
    }
    .set { pks_positive_tblouts }

  // ── Prokka annotation (pks+ bins only) ───────────────────────────────
  // Join pks+ tblouts back to their bin FASTA for Prokka input
  per_bin
    .map { meta, bin -> tuple(meta.id + "." + bin.baseName, meta, bin) }
    .join(
      pks_positive_tblouts.map { meta, tbl ->
        tuple(meta.id + "." + tbl.baseName.replaceAll(/\.tblout$/, ''), meta, tbl)
      }
    )
    .map { key, meta, bin, _meta2, tbl -> tuple(meta, bin, tbl) }
    .set { pks_bins_with_tblout }

  PROKKA(pks_bins_with_tblout.map { meta, bin, tbl -> tuple(meta, bin) }, [], [])

  // ── Genomic context extraction ─────────────────────────────────────────
  PROKKA.out.gff
    .join(pks_bins_with_tblout.map { meta, bin, tbl -> tuple(meta, tbl) })
    .map { meta, gff, tblout ->
      def sampleID = meta.id.split(/\./)[0]
      def binID    = meta.id.split(/\./, 2)[1] ?: meta.id
      tuple(sampleID, binID, gff, tblout)
    }
    .set { context_input }

  extractGenomicContext(context_input)

  // ── Collect all outputs → cohort summary table ────────────────────────
  CHECKM2_PREDICT.out.output
    .map { meta, report -> report }
    .collect()
    .set { all_checkm2 }

  GTDBTK_CLASSIFYWF.out.output
    .map { meta, report -> report }
    .collect()
    .set { all_gtdbtk }

  HMMER_HMMSEARCH.out.target_summary
    .map { meta, tbl -> tbl }
    .collect()
    .set { all_tblouts }

  extractGenomicContext.out
    .map { sampleID, binID, ctx -> ctx }
    .collect()
    .set { all_context }

  magSummaryTable(all_checkm2, all_gtdbtk, all_tblouts, all_context)

  emit:
  summary = magSummaryTable.out
}
```

- [ ] **Step 3: Verify Nextflow syntax with stub-run (see Task 9)**

This step cannot be fully validated until main.nf is wired (Task 8). Proceed to Task 8 now.

- [ ] **Step 4: Commit**

```bash
git add Modules/pks_mag.nf
git commit -m "Add pksMAG workflow to pks_mag.nf"
```

---

## Task 8: Wire MAG branch into main.nf

**Files:**
- Modify: `main.nf`

Three changes: (1) new params, (2) include, (3) MAG branch block at the end of the workflow.

- [ ] **Step 1: Add params to main.nf**

After the existing `params.kraken_db = null` line (line ~40), add:

```groovy
params.enable_mags     = false
params.gtdbtk_db       = null
params.checkm2_db      = null
params.clb_protein_hmm = "${projectDir}/ref/hmm/clb_all_protein.hmm"
```

- [ ] **Step 2: Add include after existing includes (~line 65)**

```groovy
include { pksMAG } from './Modules/pks_mag.nf'
```

- [ ] **Step 3: Add validation in the workflow block**

After the existing `pks_taxa` validation block (around line 89), add:

```groovy
if (params.enable_mags && params.input_context != "metagenome") {
    exit 1, "--enable_mags requires --input_context metagenome"
}
if (params.enable_mags && !params.gtdbtk_db) {
    exit 1, "--enable_mags requires --gtdbtk_db"
}
if (params.enable_mags && !params.checkm2_db) {
    exit 1, "--enable_mags requires --checkm2_db"
}
```

- [ ] **Step 4: Add MAG branch at the end of the workflow block**

Immediately before the final `}` closing the `workflow` block, after the `masterQCSummary` call:

```groovy
// ---------- STEP 6: MAG branch (optional) ----------
if (params.enable_mags && params.input_context == "metagenome") {
    pksMAG(MAPPED_READS)
}
```

- [ ] **Step 5: Commit**

```bash
git add main.nf
git commit -m "Wire MAG branch into main.nf (params + include + channel)"
```

---

## Task 9: Validate with stub-run

**Files:**
- None changed

A Nextflow stub-run parses all DSL2 syntax and channel wiring without executing any process. It confirms the workflow compiles and the MAG branch is reachable.

- [ ] **Step 1: Run stub-run without MAG (baseline — must still pass)**

From an interactive platinum node with `env_nf` active:

```bash
cd /tscc/lustre/restricted/alexandrov-ddn/users/l1joseph/pksProfiler/runs
mkdir -p stub_test && cd stub_test

nextflow run ../../main.nf \
  -stub-run \
  --sample ../../docs/cm017_smoke_test.csv \
  --input_context metagenome \
  --input_data_type fastq \
  --hg38_db dummy \
  --t2t_phix_db dummy \
  2>&1 | tail -20
```

Expected: workflow completes, no `pksMAG` jobs appear (enable_mags is false).

- [ ] **Step 2: Run stub-run with MAG branch enabled**

```bash
nextflow run ../../main.nf \
  -stub-run \
  --sample ../../docs/cm017_smoke_test.csv \
  --input_context metagenome \
  --input_data_type fastq \
  --hg38_db dummy \
  --t2t_phix_db dummy \
  --enable_mags true \
  --gtdbtk_db /tscc/lustre/restricted/alexandrov-ddn/users/l1joseph/databases/gtdbtk \
  --checkm2_db dummy \
  2>&1 | tail -30
```

Expected output includes process names from the MAG branch:
```
[xx/xxxxxx] process > MEGAHIT (CM017)                             [100%]
[xx/xxxxxx] process > MINIMAP2_ALIGN (CM017)                      [100%]
[xx/xxxxxx] process > METABAT2_METABAT2 (CM017)                   [100%]
[xx/xxxxxx] process > CHECKM2_PREDICT (CM017)                     [100%]
[xx/xxxxxx] process > GTDBTK_CLASSIFYWF (CM017)                   [100%]
[xx/xxxxxx] process > PRODIGAL (CM017.bin.X)                      [100%]
[xx/xxxxxx] process > HMMER_HMMSEARCH (CM017.bin.X)               [100%]
[xx/xxxxxx] process > magSummaryTable                             [100%]
```

- [ ] **Step 3: Fix any channel wiring errors**

If the stub-run fails with a channel shape mismatch, read the error, find the process in `Modules/pks_mag.nf` with the wrong input tuple, and adjust. Common issues:
- nf-core module expected `tuple val(meta), path(ref)` but got positional args → wrap reference in `tuple([id: 'ref'], file(...))`
- `flatMap` on single file (not list) → add `(bins instanceof List ? bins : [bins])`

- [ ] **Step 4: Commit any fixes**

```bash
git add Modules/pks_mag.nf
git commit -m "Fix pksMAG channel wiring from stub-run"
```

---

## Self-Review Checklist

**Spec coverage:**
- ✅ Protein HMM build → Task 1
- ✅ nf-core modules installed → Task 3
- ✅ Resource labels → Task 2
- ✅ MEGAHIT assembly → pksMAG workflow Task 7
- ✅ minimap2 → contig depth → MetaBAT2 binning → Task 7
- ✅ CheckM2 + GTDB-Tk → Task 7
- ✅ Prodigal on ALL bins → Task 7
- ✅ hmmsearch clb_all_protein.hmm (no taxonomy pre-filter) → Task 7
- ✅ pks+ bin filter → Task 7
- ✅ Prokka (pks+ only) → Task 7
- ✅ extract_genomic_context.py → Task 4
- ✅ pks_mag_summary.tsv schema (all 13 columns) → Task 5
- ✅ magSummaryTable process → Task 6
- ✅ main.nf wiring + validation guards → Task 8
- ✅ Stub-run test → Task 9
- ✅ unexpected_taxon_flag → build_mag_summary.py
- ✅ GTDB-Tk existing DB at `/tscc/.../databases/gtdbtk` → params default in Task 8

**Not in scope (deferred):** DIAMOND rescue, TaxID 543→91347, final evidence classifier.
