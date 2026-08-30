# Enterobacteriaceae Pre-filter (v2) — Design Spec

**Date:** 2026-08-29
**Branch:** upstream-main
**Status:** Approved, implementation pending

---

## Problem

`pksProfilerHMM` and `pksProfilerAlign` currently receive all ~40M host-depleted reads per sample. HMM runtime is ~11h per sample as a result. In the human gut, pks island signal is exclusively from Enterobacteriaceae (*E. coli*, *Klebsiella*, *Citrobacter*, *Enterobacter*). Filtering reads to Enterobacteriaceae before profiling reduces the search space to ~4k–200k reads with no loss of true pks+ signal.

---

## Architecture

### Current flow
```
filterReads → mapReads → MAPPED_READS → pksProfilerAlign
                                      → pksProfilerHMM
                                      → extractPksIslandReads (via ALIGN_OUT)
```

### New flow
```
filterReads → mapReads → MAPPED_READS (QC only)
                               │
                               ▼
                    filterEnterobacteriaceae (new)
                               │
                               ▼
                        PROFILING_READS ──→ pksProfilerAlign
                                        ──→ pksProfilerHMM
                                        ──→ extractPksIslandReads (via ALIGN_OUT)
```

`MAPPED_READS` is retained solely for its existing QC emit (`reads_after_hg38`, `reads_after_t2t_phix`). `PROFILING_READS` is the filtered channel passed to all profiling steps.

---

## New Parameter

```groovy
params.entero_filter = true
```

- Default `true` — intended as the standard way to run the pipeline in v2.
- If `params.kraken_db` is null, the filter is silently skipped with a loud `log.warn`. Pipeline continues using `MAPPED_READS` as `PROFILING_READS`.
- No new `--kraken_db` validation beyond what already exists for `--pks_taxa`.

---

## New Process: `filterEnterobacteriaceae`

**File:** `Modules/pks_taxa.nf` (appended — keeps all KrakenUniq logic together)

**Resource label:** `entero_filter` (new label in `conf/base.config`)
- cpus: 8
- memory: 256 GB (KrakenUniq DB loading)
- time: 24h

**Conda:** `krakenuniq_bracken_env` (already has `krakenuniq` + `krakentools`)

**scratch:** `true`

### Inputs
```
tuple val(sampleID), path(fastq_gz)   ← MAPPED_READS
```

### Outputs
```
tuple val(sampleID), path("${sampleID}.entero.fastq.gz"),  emit: reads
tuple val(sampleID), path("${sampleID}.entero.qc.tsv"),    emit: qc
```

### Script logic
1. Decompress input FASTQ.
2. **Empty input guard:** if FASTQ is empty, write empty gzip output + zero-count QC TSV and exit 0.
3. Run `krakenuniq --db --threads --report-file --output` on all reads.
4. Run `extract_kraken_reads.py -t 543 --include-children` using report + output files to extract Enterobacteriaceae reads.
5. gzip filtered FASTQ → `${sampleID}.entero.fastq.gz`.
6. Count reads before and after filter → write QC TSV:
   ```
   Sample  Metric                      Value
   CM017   reads_after_entero_filter   4231
   ```
7. **Empty output guard:** if zero reads survive, write empty gzip + zero-count QC and exit 0 (same pattern as existing modules).

---

## `main.nf` Changes

### 1. Import
Add `filterEnterobacteriaceae` to the existing `pks_taxa.nf` include:
```groovy
include { filterEnterobacteriaceae; extractPksIslandReads; Bracken; ... } from './Modules/pks_taxa.nf'
```

### 2. Parameter declaration
Add to params section:
```groovy
params.entero_filter = true
```

### 3. Channel wiring (after `mapReads`)
```groovy
// STEP 1c: Enterobacteriaceae pre-filter
def PROFILING_READS
if (params.entero_filter && params.kraken_db) {
    FILTER_ENTERO_OUT = filterEnterobacteriaceae(MAPPED_READS)
    PROFILING_READS = FILTER_ENTERO_OUT.reads
    QC_FRAGMENTS = QC_FRAGMENTS.mix(
        FILTER_ENTERO_OUT.qc.map { _id, qc -> qc }
    )
} else {
    if (params.entero_filter && !params.kraken_db) {
        log.warn "entero_filter is enabled but --kraken_db not provided. " +
                 "Profiling will run on all host-depleted reads (~40M). " +
                 "Provide --kraken_db to enable the Enterobacteriaceae pre-filter."
    }
    PROFILING_READS = MAPPED_READS
}
```

### 4. Swap profiling inputs
- `pksProfilerAlign(MAPPED_READS)` → `pksProfilerAlign(PROFILING_READS)`
- `pksProfilerHMM(MAPPED_READS)` → `pksProfilerHMM(PROFILING_READS)`
- `extractPksIslandReads` receives BAM from `ALIGN_OUT` — unchanged.

### 5. Lift HMM restriction
Remove:
```groovy
if (params.pks_taxa && params.profiling_method == "hmm") {
    exit 1, "--pks_taxa requires alignment profiling. ..."
}
```
No longer valid — HMM now receives pre-filtered reads regardless.

---

## `conf/base.config` Changes

Add new label:
```groovy
withLabel:entero_filter {
    cpus   = 8
    memory = { 256.GB * task.attempt }
    time   = { 24.h  * task.attempt }
}
```

---

## Files Changed

| File | Change |
|------|--------|
| `Modules/pks_taxa.nf` | Add `filterEnterobacteriaceae` process |
| `main.nf` | Import, param, channel wiring, swap profiling inputs, lift HMM restriction |
| `conf/base.config` | Add `entero_filter` resource label |

No changes to: `pksProfiler_hmm.nf`, `pksProfiler_align.nf`, `conf/tscc.config`.

---

## Testing

1. **Stub run** — `nextflow run main.nf -stub-run -profile tscc --kraken_db <db>` — confirms parsing and channel wiring, no Groovy/DSL errors.
2. **Single-sample smoke test** — one Shiba 2026 sample with `--entero_filter true --kraken_db <db>`: verify `entero.fastq.gz` produced, `reads_after_entero_filter` in QC summary, profiling steps receive filtered FASTQ.
3. **No-DB fallback** — run with `--entero_filter true` but no `--kraken_db`: confirm WARN in `.nextflow.log`, pipeline completes on full `MAPPED_READS`.
4. **Correctness** — CM017 (known pks+, 59k reads): alignment read counts identical with filtered vs unfiltered input.

---

## Notes

- No PR to upstream for now — for local testing only.
- Ammal is aware of this design; coordinate before merging to `upstream-main`.
- KrakenUniq chosen over Kraken2 for consistency with existing `pks_taxa` module. Kraken2 is a valid future upgrade for better recall on divergent strains.
