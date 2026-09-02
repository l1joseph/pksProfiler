# pksProfiler v3: MAG Module Design Spec

**Goal:** Add an optional metagenome-assembled genome (MAG) branch to pksProfiler that assembles bins from clean non-host reads, detects pks+ MAGs using colibactin protein HMMs, assigns taxonomy via GTDB-Tk, and reports genomic context (integrases, transposases, flanking synteny) for each pks+ bin.

**Architecture:** Two parallel branches downstream of host depletion — a read-level branch (existing) and a new MAG branch (`--enable_mags true`, metagenome mode only) — both feeding into a unified evidence summary.

**Tech Stack:** nf-core DSL2 modules (MEGAHIT, minimap2, MetaBAT2, CheckM2, GTDB-Tk, Prodigal, HMMER, Prokka), MAFFT + HMMER for offline protein HMM build, DIAMOND for offline db build, Python for summary and genomic context scripts.

**Implementation scope:** `Modules/pks_mag.nf` and MAG summary table process only. The full v3 read-level upgrades (DIAMOND rescue, TaxID 91347) are documented below but deferred.

---

## Full v3 Pipeline (complete vision)

### Preprocessing — unchanged

```
BAM / CRAM / FASTQ
  → extractReads (samtools)
  → filterReads (fastp)
  → mapReads (minimap2: hg38 → T2T → phiX)
  → MAPPED_READS  (all clean non-host reads)
```

### Read-level branch — v3 upgrade (deferred except TaxID fix)

Upgrade from current `filterEnterobacteriaceae` (TaxID 543) to full v3 read rescue:

```
MAPPED_READS
  → KrakenUniq (TaxID 91347 Enterobacterales)
      ├─ classified → primary candidates + both mates
      └─ unclassified → DIAMOND blastx (curated ClbA–ClbS protein db)
                             └─ hits → rescue candidates + both mates
  → merge + deduplicate → PROFILING_READS
  → Bowtie2 → IHE3034 single-strain index (reference-like evidence)
  → nhmmscan DNA HMMs — clb_all_dna.hmm (divergent evidence)
  → Bracken (genus/species abundance)
```

**Design notes:**
- TaxID 91347 (Enterobacterales, order) replaces 543 (Enterobacteriaceae, family) — broader to catch Yersinia, Serratia, etc.
- Bowtie2 stays single-strain (IHE3034): nhmmscan and DIAMOND handle divergent sequences
- DIAMOND db: same curated ClbA–ClbS protein set used for protein HMM build
- Existing `pksProfilerAlign` + `pksProfilerHMM` + `plotting.nf` + all scripts unchanged

**Reuse from current pipeline:**

| Component | Status |
|-----------|--------|
| `extract_reads.nf`, `filter_reads.nf`, `map_reads.nf` | Unchanged |
| `pksProfiler_align.nf`, `pksProfiler_hmm.nf` | Unchanged |
| `plotting.nf`, all scripts | Unchanged |
| `ref/hmm/clb_all_dna.hmm`, IHE3034 index, GFF | Unchanged |
| `pks_taxa.nf` — filterEnterobacteriaceae | TaxID 543 → 91347 (1 line, deferred) |
| `pks_taxa.nf` — Bracken, process_bracken, combineClbTaxonomySupport | Unchanged |

---

### MAG branch — implementation scope

Triggered by `params.enable_mags = true` AND `params.input_context == "metagenome"`.

Input: `MAPPED_READS` (all clean non-host reads, NOT entero-filtered — assembly benefits from full community).

```
MAPPED_READS
  │
  ├─ MEGAHIT assembly → contigs.fa
  │
  ├─ minimap2 align reads → contigs.fa → BAM
  ├─ samtools sort/index BAM
  ├─ jgi_summarize_bam_contig_depths → depth.txt
  │
  ├─ MetaBAT2 binning (depth.txt + contigs.fa) → bin_*.fa + unbinned.fa
  │
  ├─ CheckM2 predict (all bins) → quality_report.tsv
  │     [completeness, contamination, genome_size]
  │
  ├─ GTDB-Tk classifywf (all bins) → gtdbtk.bac120.summary.tsv
  │     [DB: params.gtdbtk_db, r220, 107 GB, already on disk]
  │
  ├─ Prodigal -p meta (ALL bins + unbinned contigs, no taxonomy pre-filter)
  │     → per-bin .faa (predicted proteins)
  │
  ├─ hmmer/hmmsearch vs clb_all_protein.hmm (ALL bins)
  │     [e-value ≤ 1e-5; catches clbS HGT in non-Enterobacteriaceae]
  │     → per-bin tblout
  │
  ├─ [filter] pks+ bins: ≥1 clb gene hit
  │
  ├─ Prokka (pks+ bins only) → .gff annotation
  │
  ├─ extract_genomic_context.py (pks+ bins)
  │     ├─ extract ±50 kb flanking region around each clb hit
  │     ├─ hmmsearch Pfam HMMs: integrases (PF00589, PF13102), transposases (PF01526)
  │     └─ report flanking gene names (synteny)
  │
  └─ pks_mag_summary process → pks_mag_summary.tsv
```

**Key design decisions:**
- `MAPPED_READS` (not `PROFILING_READS`) as assembly input — full community coverage improves binning
- nf-core modules provide versioned, tested tool wrappers; channel adapter converts `tuple val(sampleID), path(reads)` → `tuple val(meta), path(reads)` where `meta = [id: sampleID]`
- Prodigal and hmmsearch run on ALL bins regardless of GTDB-Tk taxonomy — clbS homologs in Bacteroides/Firmicutes (HGT) would be missed if filtered
- `unexpected_taxon_flag`: set true for any pks+ bin assigned outside Enterobacterales by GTDB-Tk
- Unbinned contigs searched alongside bins — Ammal's v3 diagram explicitly marks this

---

## nf-core Modules to Install

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

Installed to `modules/nf-core/<tool>/main.nf`. Each carries its own versioned conda/container spec — no new `conda_envs/` YML required.

---

## Prerequisite: Protein HMM Build (offline, one-time)

Build `ref/hmm/clb_all_protein.hmm` before the MAG module can run.

**Sources:**
- Primary: MIBiG BGC0000943 (colibactin, E. coli IHE3034) — all 19 clbA–clbS proteins
- Additional strains for diversity (5–10 total):
  - E. coli Nissle 1917 (B2 phylogroup, probiotic)
  - E. coli CFT073 (B2, UPEC)
  - E. coli 536 (B2, UPEC)
  - E. coli APEC O1 (avian pathogenic)
  - Klebsiella pneumoniae pks+ strain (cross-genus)
  - Enterobacter cloacae (if annotated clb available)
- **clbS specifically:** also include any non-Enterobacteriaceae clbS-like homologs from literature (resistance gene HGT candidates)

**Build process (per gene):**
```bash
# 1. Fetch proteins for each clb gene from NCBI across strains
esearch -db protein -query "clbA[Gene Name] Escherichia coli" | efetch -format fasta > clbA.fasta
# repeat for clbB..clbS

# 2. Multiple sequence alignment
mafft --auto clbA.fasta > clbA.msa

# 3. Build per-gene HMM
hmmbuild clbA.hmm clbA.msa

# 4. Concatenate all 19 HMMs
cat clb*.hmm > clb_all_protein.hmm

# 5. Press for hmmsearch
hmmpress clb_all_protein.hmm
```

Committed to `ref/hmm/clb_all_protein.hmm` + `.h3{f,i,m,p}`.

**DIAMOND db (for deferred DIAMOND rescue step):**
```bash
cat clb*.fasta > clb_all_proteins.fasta
diamond makedb --in clb_all_proteins.fasta --db ref/diamond/clb_proteins
```

---

## New Parameters

```groovy
params.enable_mags     = false
params.gtdbtk_db       = "/tscc/lustre/restricted/alexandrov-ddn/users/l1joseph/databases/gtdbtk"
params.clb_protein_hmm = "${projectDir}/ref/hmm/clb_all_protein.hmm"
params.diamond_clb_db  = "${projectDir}/ref/diamond/clb_proteins.dmnd"  // deferred
// MAG outputs land under params.outdir/pks_summary/mag/ via publishDir in each process
```

---

## Outputs

```
results/pks_summary/mag/
├─ pks_mag_summary.tsv         # cohort-level table (one row per pks+ bin)
├─ <sample>/<bin_id>.faa       # predicted proteins (pks+ bins)
├─ <sample>/<bin_id>.gff       # Prokka annotation (pks+ bins)
└─ <sample>/<bin_id>.context/  # genomic context files
```

### `pks_mag_summary.tsv` schema

| Column | Description |
|--------|-------------|
| `sample` | Sample ID |
| `bin_id` | MetaBAT2 bin identifier |
| `taxonomy` | GTDB-Tk classification (full lineage) |
| `completeness` | CheckM2 completeness % |
| `contamination` | CheckM2 contamination % |
| `clb_genes_detected` | Count of distinct clb genes with hits |
| `clb_genes` | Comma-separated list (e.g. clbA,clbB,clbN) |
| `best_evalue` | Best hmmsearch e-value across all hits |
| `has_integrase` | Boolean — Pfam integrase HMM hit in ±50kb |
| `has_transposase` | Boolean — Pfam transposase HMM hit in ±50kb |
| `flanking_genes` | Semicolon-separated flanking gene names (Prokka) |
| `unexpected_taxon_flag` | True if GTDB-Tk assigns outside Enterobacterales |

---

## Resource Labels (conf/base.config additions)

```groovy
withLabel: mag_assembly {
    cpus   = 16
    memory = { 64.GB * task.attempt }
    time   = { 12.h  * task.attempt }
}
withLabel: mag_binning {
    cpus   = 8
    memory = { 32.GB * task.attempt }
    time   = { 4.h   * task.attempt }
}
withLabel: mag_gtdbtk {
    cpus   = 16
    memory = { 200.GB * task.attempt }
    time   = { 24.h  * task.attempt }
}
withLabel: mag_hmm {
    cpus   = 8
    memory = { 16.GB * task.attempt }
    time   = { 4.h   * task.attempt }
}
```

---

## Files Changed / Created

| File | Action |
|------|--------|
| `Modules/pks_mag.nf` | **Create** — wrapper workflow + custom processes |
| `modules/nf-core/*/main.nf` | **Install** via nf-core CLI |
| `scripts/extract_genomic_context.py` | **Create** |
| `scripts/build_mag_summary.py` | **Create** |
| `ref/hmm/clb_all_protein.hmm` | **Build offline, commit** |
| `ref/diamond/clb_proteins.dmnd` | **Build offline, commit** (deferred) |
| `conf/base.config` | **Add** 4 new resource labels |
| `main.nf` | **Add** `enable_mags` param + import + MAG branch wiring |

**Not changed:** `pksProfiler_align.nf`, `pksProfiler_hmm.nf`, `plotting.nf`, all existing scripts, all existing ref files.

---

## Deferred (future work)

- TaxID 543 → 91347 in `filterEnterobacteriaceae` (1-line change, low risk, defer to coordinate with Ammal)
- DIAMOND blastx rescue step (requires curated protein db + new process in pks_taxa.nf)
- Final evidence classifier combining read-level + MAG evidence
- Genomic context: full synteny graph beyond ±50kb flanking
