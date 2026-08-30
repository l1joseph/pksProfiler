# Enterobacteriaceae Pre-filter Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a `filterEnterobacteriaceae` process that runs KrakenUniq on all host-depleted reads right after host depletion, filters to Enterobacteriaceae (TaxID 543), and passes the filtered FASTQ as `PROFILING_READS` to all downstream profiling modules — replacing the current pattern where ~40M raw reads flow into every profiling step.

**Architecture:** New process `filterEnterobacteriaceae` in `Modules/pks_taxa.nf` emits a filtered FASTQ and a QC TSV. `main.nf` sets a `PROFILING_READS` channel to either the filtered output (when `kraken_db` is provided) or `MAPPED_READS` with a warning (when not). All three profiling steps (`pksProfilerAlign`, `pksProfilerHMM`, `extractPksIslandReads`) receive `PROFILING_READS` instead of `MAPPED_READS`.

**Tech Stack:** Nextflow DSL2, KrakenUniq, KrakenTools (`extract_kraken_reads.py`), bash, gzip/seqtk

---

## File Map

| File | Change |
|------|--------|
| `Modules/pks_taxa.nf` | Append new `filterEnterobacteriaceae` process |
| `main.nf` | Add param, import, channel wiring, swap profiling inputs, lift HMM restriction |
| `conf/base.config` | Add `entero_filter` resource label |

---

### Task 1: Add `entero_filter` resource label to `conf/base.config`

**Files:**
- Modify: `conf/base.config`

- [ ] **Step 1: Open `conf/base.config` and locate the last `withLabel` block**

The file currently ends with:
```groovy
    withLabel:pks_hmm {
        cpus   = 8
        memory = { 128.GB * task.attempt }
        time   = { 30.h * task.attempt }
    }
}
```

- [ ] **Step 2: Append the new label inside the `process { }` block, before the closing `}`**

Final state of the last two labels:
```groovy
    withLabel:pks_hmm {
        cpus   = 8
        memory = { 128.GB * task.attempt }
        time   = { 30.h * task.attempt }
    }
    withLabel:entero_filter {
        cpus   = 8
        memory = { 256.GB * task.attempt }
        time   = { 24.h  * task.attempt }
    }
}
```

- [ ] **Step 3: Verify syntax**

```bash
nextflow config -profile tscc 2>&1 | head -5
```
Expected: no `Error` lines, config dumps cleanly.

- [ ] **Step 4: Commit**

```bash
git add conf/base.config
git commit -m "Add entero_filter resource label to base.config"
```

---

### Task 2: Add `filterEnterobacteriaceae` process to `Modules/pks_taxa.nf`

**Files:**
- Modify: `Modules/pks_taxa.nf` (append after line 336)

- [ ] **Step 1: Append the new process at the end of `Modules/pks_taxa.nf`**

Add the following after the last line of the file (after `combineClbTaxonomySupport`):

```groovy


process filterEnterobacteriaceae {
  label 'entero_filter'
  scratch true
  publishDir "${params.pks_dir}", mode: 'copy'
  conda "${params.krakenuniq_bracken_env}"

  input:
  tuple val(sampleID), path(fastq_gz)

  output:
  tuple val(sampleID), path("${sampleID}.entero.fastq.gz"), emit: reads
  tuple val(sampleID), path("${sampleID}.entero.qc.tsv"),   emit: qc

  script:
  """
  set -euo pipefail

  REPORT="${sampleID}.krakenuniq.entero.report.txt"
  OUTPUT="${sampleID}.krakenuniq.entero.output.txt"
  FILTERED="${sampleID}.entero.fastq"
  QC="${sampleID}.entero.qc.tsv"

  zcat "${fastq_gz}" > "${sampleID}.all_reads.fastq"

  READS_IN=\$(awk 'END { print int(NR / 4) }' "${sampleID}.all_reads.fastq")

  printf "Sample\tMetric\tValue\n" > "\$QC"

  # Empty input guard: skip KrakenUniq and write empty outputs
  if [[ "\$READS_IN" -eq 0 ]]; then
    : | gzip -c > "${sampleID}.entero.fastq.gz"
    printf "%s\treads_after_entero_filter\t0\n" "${sampleID}" >> "\$QC"
    exit 0
  fi

  krakenuniq \\
    --db "${params.kraken_db}" \\
    --threads "${task.cpus}" \\
    --report-file "\$REPORT" \\
    --output "\$OUTPUT" \\
    "${sampleID}.all_reads.fastq"

  # Extract reads classified under Enterobacteriaceae (TaxID 543) and all descendants
  extract_kraken_reads.py \\
    -k "\$OUTPUT" \\
    -r "\$REPORT" \\
    -s "${sampleID}.all_reads.fastq" \\
    -t 543 \\
    --include-children \\
    -o "\$FILTERED"

  # Empty output guard
  if [[ ! -s "\$FILTERED" ]]; then
    : | gzip -c > "${sampleID}.entero.fastq.gz"
    printf "%s\treads_after_entero_filter\t0\n" "${sampleID}" >> "\$QC"
    exit 0
  fi

  gzip -c "\$FILTERED" > "${sampleID}.entero.fastq.gz"

  READS_OUT=\$(awk 'END { print int(NR / 4) }' "\$FILTERED")
  printf "%s\treads_after_entero_filter\t%s\n" "${sampleID}" "\$READS_OUT" >> "\$QC"
  """
}
```

- [ ] **Step 2: Verify the process is syntactically valid**

```bash
nextflow run main.nf -stub-run -profile tscc \
  --sample docs/shiba2026_samples.csv \
  --input_data_type fastq \
  --kraken_db /tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/krakenUniq_8_8_2023 \
  2>&1 | head -30
```
Expected: process list includes `filterEnterobacteriaceae`, no `No such variable` or `Unexpected token` errors.

- [ ] **Step 3: Commit**

```bash
git add Modules/pks_taxa.nf
git commit -m "Add filterEnterobacteriaceae process to pks_taxa.nf"
```

---

### Task 3: Wire `filterEnterobacteriaceae` into `main.nf`

**Files:**
- Modify: `main.nf`

#### Step 3a — Add `params.entero_filter`

- [ ] **Step 1: Add parameter declaration after line 12 (`params.hmm_chunking = false`)**

Current block (lines 10–13):
```groovy
params.profiling_method = "bowtie2" // bowtie2 | hmm | both
params.hmm_evalue       = 1e-10
params.hmm_chunking = false
params.hmm_model        = "${projectDir}/ref/hmm/clb_all_dna.hmm"
```

Replace with:
```groovy
params.profiling_method  = "bowtie2" // bowtie2 | hmm | both
params.hmm_evalue        = 1e-10
params.hmm_chunking      = false
params.hmm_model         = "${projectDir}/ref/hmm/clb_all_dna.hmm"
params.entero_filter     = true      // filter host-depleted reads to Enterobacteriaceae before profiling
```

#### Step 3b — Add `filterEnterobacteriaceae` to the import

- [ ] **Step 2: Update the `pks_taxa.nf` include on line 62**

Current:
```groovy
include { extractPksIslandReads; Bracken; process_bracken as combinePKSTaxa; combineClbTaxonomySupport } from './Modules/pks_taxa.nf'
```

Replace with:
```groovy
include { filterEnterobacteriaceae; extractPksIslandReads; Bracken; process_bracken as combinePKSTaxa; combineClbTaxonomySupport } from './Modules/pks_taxa.nf'
```

#### Step 3c — Insert channel wiring after `mapReads`

- [ ] **Step 3: Insert the `PROFILING_READS` block after the `MAP_OUT` section**

Current (lines 216–223):
```groovy
	MAP_OUT = mapReads(FILTERED_UNMAPPED_READS)

	MAP_OUT.reads
	    .set { MAPPED_READS }

    QC_FRAGMENTS = QC_FRAGMENTS.mix(
        MAP_OUT.qc.map { _sampleID, qc_file -> qc_file }
    )
```

Replace with:
```groovy
	MAP_OUT = mapReads(FILTERED_UNMAPPED_READS)

	MAP_OUT.reads
	    .set { MAPPED_READS }

    QC_FRAGMENTS = QC_FRAGMENTS.mix(
        MAP_OUT.qc.map { _sampleID, qc_file -> qc_file }
    )

    // ---------- STEP 1c: Enterobacteriaceae pre-filter ----------
    def PROFILING_READS
    if (params.entero_filter && params.kraken_db) {
        FILTER_ENTERO_OUT = filterEnterobacteriaceae(MAPPED_READS)
        PROFILING_READS = FILTER_ENTERO_OUT.reads
        QC_FRAGMENTS = QC_FRAGMENTS.mix(
            FILTER_ENTERO_OUT.qc.map { _sampleID, qc_file -> qc_file }
        )
    } else {
        if (params.entero_filter && !params.kraken_db) {
            log.warn "[pksProfiler] entero_filter is enabled but --kraken_db not provided. " +
                     "Profiling will run on all host-depleted reads (~40M). " +
                     "Provide --kraken_db to enable the Enterobacteriaceae pre-filter."
        }
        PROFILING_READS = MAPPED_READS
    }
```

#### Step 3d — Swap `MAPPED_READS` → `PROFILING_READS` in profiling calls

- [ ] **Step 4: Update `pksProfilerAlign` call**

Current:
```groovy
        ALIGN_OUT = pksProfilerAlign(MAPPED_READS)
```
Replace with:
```groovy
        ALIGN_OUT = pksProfilerAlign(PROFILING_READS)
```

- [ ] **Step 5: Update `pksProfilerHMM` call**

Current:
```groovy
        HMM_OUT = pksProfilerHMM(MAPPED_READS)
```
Replace with:
```groovy
        HMM_OUT = pksProfilerHMM(PROFILING_READS)
```

#### Step 3e — Lift the HMM + pks_taxa restriction

- [ ] **Step 6: Remove the exit block that prevents running `pks_taxa` with HMM**

Current (lines 232–234):
```groovy
    if (params.pks_taxa && params.profiling_method == "hmm") {
        exit 1, "--pks_taxa requires alignment profiling. Use --profiling_method bowtie2 or both."
    }
```
Delete these three lines entirely.

- [ ] **Step 7: Verify stub run with all three modes**

```bash
nextflow run main.nf -stub-run -profile tscc \
  --sample docs/shiba2026_samples.csv \
  --input_data_type fastq \
  --profiling_method both \
  --pks_taxa true \
  --kraken_db /tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/krakenUniq_8_8_2023 \
  --bracken_read_length 150 \
  2>&1 | head -40
```
Expected: `filterEnterobacteriaceae` appears in the process list, no exit errors, `pks_taxa` + HMM combination no longer blocked.

- [ ] **Step 8: Commit**

```bash
git add main.nf
git commit -m "Wire Enterobacteriaceae pre-filter into main.nf"
```

---

### Task 4: Smoke test on a single known pks+ sample

**Files:**
- Read: `runs/shiba2026_v2/work/<hash>/.command.sh` (after run)

- [ ] **Step 1: Launch single-sample test with filter enabled**

From an interactive platinum node in `runs/shiba2026_v2/`:
```bash
nextflow run ../../main.nf \
  -profile tscc \
  --sample <(head -2 /tscc/lustre/restricted/alexandrov-ddn/users/l1joseph/pksProfiler/docs/shiba2026_samples.csv | grep CM017 || echo "patient,fastq1,fastq2
CM017,$(grep CM017 /tscc/lustre/restricted/alexandrov-ddn/users/l1joseph/pksProfiler/docs/shiba2026_samples.csv | cut -d, -f2),$(grep CM017 /tscc/lustre/restricted/alexandrov-ddn/users/l1joseph/pksProfiler/docs/shiba2026_samples.csv | cut -d, -f3)") \
  --input_data_type fastq \
  --profiling_method bowtie2 \
  --kraken_db /tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/krakenUniq_8_8_2023 \
  --outdir ../../RESULTS/entero_filter_test \
  2>&1 | tee test_entero_filter.log
```

Or simpler — create a one-sample CSV first:
```bash
head -1 docs/shiba2026_samples.csv > /tmp/cm017.csv
grep "^CM017," docs/shiba2026_samples.csv >> /tmp/cm017.csv
```
Then:
```bash
nextflow run main.nf \
  -profile tscc \
  --sample /tmp/cm017.csv \
  --input_data_type fastq \
  --profiling_method bowtie2 \
  --kraken_db /tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/krakenUniq_8_8_2023 \
  --outdir RESULTS/entero_filter_test \
  2>&1 | tee test_entero_filter.log
```

- [ ] **Step 2: Verify `entero.fastq.gz` was produced with non-zero reads**

```bash
zcat RESULTS/entero_filter_test/pks_per_sample/CM017.entero.fastq.gz | awk 'END { print int(NR/4), "reads" }'
```
Expected: non-zero read count, substantially less than ~40M (should be hundreds to low thousands for a pks+ sample).

- [ ] **Step 3: Verify `reads_after_entero_filter` appears in QC summary**

```bash
grep "entero_filter" RESULTS/entero_filter_test/pks_summary/qc/cohort_qc_summary.tsv
```
Expected: one row for CM017 with `reads_after_entero_filter` metric.

- [ ] **Step 4: Verify `pksProfilerAlign` received the filtered FASTQ**

Find the `pksProfilerAlign` work dir:
```bash
grep "pksProfilerAlign" test_entero_filter.log | grep -oP '\[\w+/\w+\]' | head -1
```
Then check its command:
```bash
cat runs/shiba2026_v2/work/<hash>/.command.sh | grep "entero"
```
Expected: input FASTQ path ends in `.entero.fastq.gz`, not `.UNMAPPED.FASTP.FILTERED.hg38.t2t.fastq.gz`.

- [ ] **Step 5: Correctness check — alignment read counts unchanged**

CM017 alignment result with filter should match the cached result from the full run (59,555 total reads on pks island):
```bash
grep "CM017" RESULTS/entero_filter_test/pks_summary/gene_counts/pks.gene.counts.align.txt | \
  awk '{sum=0; for(i=2;i<=NF;i++) sum+=$i; print sum, "total reads"}'
```
Expected: ~59,555 (within a few reads — all pks+ reads are *E. coli* and survive the filter).

- [ ] **Step 6: Commit**

```bash
git add test_entero_filter.log   # optional, skip if not wanted in repo
git commit -m "Verified Enterobacteriaceae pre-filter smoke test on CM017"
```

---

### Task 5: No-DB fallback test

**Files:**
- Read: `.nextflow.log` after run

- [ ] **Step 1: Run without `--kraken_db`**

```bash
nextflow run main.nf \
  -profile tscc \
  --sample /tmp/cm017.csv \
  --input_data_type fastq \
  --profiling_method bowtie2 \
  --outdir RESULTS/entero_filter_nodb_test \
  2>&1 | tee test_nodb.log
```

- [ ] **Step 2: Confirm warning in log**

```bash
grep "entero_filter is enabled" test_nodb.log
```
Expected:
```
WARN: [pksProfiler] entero_filter is enabled but --kraken_db not provided. Profiling will run on all host-depleted reads (~40M). Provide --kraken_db to enable the Enterobacteriaceae pre-filter.
```

- [ ] **Step 3: Confirm `filterEnterobacteriaceae` was NOT submitted as a SLURM job**

```bash
grep "filterEnterobacteriaceae" test_nodb.log
```
Expected: no lines (process was skipped entirely).

- [ ] **Step 4: Confirm pipeline completed successfully on full reads**

```bash
grep "Succeeded\|Completed" test_nodb.log | tail -3
```
Expected: `Succeeded : N` with no failures.

- [ ] **Step 5: Commit**

```bash
git commit --allow-empty -m "Validated no-DB fallback path for entero_filter"
```
