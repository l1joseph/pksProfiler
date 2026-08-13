nextflow.enable.dsl = 2

process pksProfiler_hmm {
    label 'pksProfilerHMM'
    scratch true
    publishDir "${params.pks_dir}", mode: 'copy'
    conda "${params.pks_hmm_env}"

    input:
    tuple val(sampleID), path(r1), path(r2)

    output:
    tuple val(sampleID),
          path("${sampleID}.nhmmscan.tblout"),
          path("${sampleID}.nhmmscan.log"),
          path("${sampleID}.hmm_counts.tsv"),
          path("${sampleID}.hits.read_ids.txt"),
          path("${sampleID}.hmm.filtered.fa")

    script:
    def tblout      = "${sampleID}.nhmmscan.tblout"
    def logfile     = "${sampleID}.nhmmscan.log"
    def counts_tsv  = "${sampleID}.hmm_counts.tsv"
    def read_ids    = "${sampleID}.hits.read_ids.txt"
    def filtered_fa = "${sampleID}.hmm.filtered.fa"

    """
    set -euo pipefail

    # Reuse outputs if they already exist in publishDir
    if [[ -s "${params.pks_dir}/${tblout}" && -s "${params.pks_dir}/${logfile}" && -s "${params.pks_dir}/${counts_tsv}" && -s "${params.pks_dir}/${read_ids}" && -f "${params.pks_dir}/${filtered_fa}" ]]; then
        echo "Skipping pksProfiler_hmm: Found existing outputs in ${params.pks_dir}"
        for f in "${tblout}" "${logfile}" "${counts_tsv}" "${read_ids}" "${filtered_fa}"; do
            [[ -e "\$f" ]] || (ln -s "${params.pks_dir}/\$f" . 2>/dev/null || cp "${params.pks_dir}/\$f" .)
        done
        exit 0
    fi

    if ! gzip -t "${r1}" >/dev/null 2>&1 || ! gzip -t "${r2}" >/dev/null 2>&1; then
	  echo "WARNING: Corrupt gzip input for ${sampleID}; writing empty outputs and skipping." >&2
	  : > "${tblout}"
	  : > "${logfile}"
	  : > "${counts_tsv}"
	  : > "${read_ids}"
	  : > "${filtered_fa}"
	  exit 0
	fi

    # Merge gzipped FASTQs safely (concatenated gzip can be finicky depending on tools)
    zcat "${r1}" "${r2}" | gzip -c > "${sampleID}.merged.fastq.gz"

	ls -lh "${sampleID}.merged.fastq.gz"

    # FASTQ -> FASTA
    seqtk seq -a "${sampleID}.merged.fastq.gz" > "${sampleID}.merged.fa"

	ls -lh "${sampleID}.merged.fa"

    # Split FASTA into chunks matching the CPU allocation; each chunk runs nhmmscan --cpu 1
    N_CHUNKS=${task.cpus}
    mkdir -p "${sampleID}.fa_chunks"
    seqkit split2 -p \${N_CHUNKS} \\
        -O "${sampleID}.fa_chunks" \\
        "${sampleID}.merged.fa"

    echo "Running \${N_CHUNKS} parallel nhmmscan instances on sample: ${sampleID}"

    pids=()
    for i in \$(seq 1 \${N_CHUNKS}); do
        i_pad=\$(printf '%03d' \$i)
        chunk_fa="${sampleID}.fa_chunks/${sampleID}.merged.part_\${i_pad}.fa"
        [[ -f "\${chunk_fa}" ]] || continue
        nhmmscan --cpu 1 \\
            --tblout "${sampleID}.chunk_\${i_pad}.tblout" \\
            "${params.hmm_model}" \\
            "\${chunk_fa}" \\
            > "${sampleID}.chunk_\${i_pad}.log" 2>&1 &
        pids+=(\$!)
    done

    failed=0
    for pid in "\${pids[@]}"; do
        wait "\${pid}" || { echo "ERROR: nhmmscan chunk (pid \${pid}) failed" >&2; failed=1; }
    done
    [[ "\${failed}" -eq 0 ]] || exit 1

    # Merge tblouts: 2-line header from first chunk, data rows from all, footer from last
    first=\$(ls "${sampleID}.chunk_"*.tblout | sort | head -1)
    head -2 "\${first}" > "${tblout}"
    for f in \$(ls "${sampleID}.chunk_"*.tblout | sort); do
        grep -v '^#' "\${f}" >> "${tblout}" || true
    done
    grep '^#' "\$(ls "${sampleID}.chunk_"*.tblout | sort | tail -1)" | tail -n +3 >> "${tblout}"

    # Merge logs
    : > "${logfile}"
    for f in \$(ls "${sampleID}.chunk_"*.log | sort); do
        echo "=== \${f} ===" >> "${logfile}"
        cat "\${f}" >> "${logfile}"
    done

    # Fail early with a helpful message if tblout is empty
    if [[ ! -s "${tblout}" ]]; then
        echo "ERROR: nhmmscan did not produce a non-empty ${tblout}. See ${logfile}." >&2
        exit 1
    fi

    E="${params.hmm_evalue}"

    # Counts: best hit per query under E-value threshold
    awk -v e="\$E" '
      \$0 ~ /^#/ { next }
      NF < 14 { next }
      {
        t=\$1; q=\$3; eval=\$13+0; score=\$14+0;
        if (eval <= e) {
          if (!(q in bestE) || eval < bestE[q] || (eval == bestE[q] && score > bestScore[q])) {
            bestE[q]=eval; bestScore[q]=score; bestT[q]=t;
          }
        }
      }
      END {
        for (q in bestT) cnt[bestT[q]]++;
        for (g in cnt) printf "%s\\t%d\\n", g, cnt[g];
      }
    ' "${tblout}" | sort -k1,1 > "${counts_tsv}.tmp"

    printf "Gene\\tCount\\n" > "${counts_tsv}"
    cat "${counts_tsv}.tmp" >> "${counts_tsv}"
    rm -f "${counts_tsv}.tmp"

    # Read IDs that pass threshold
    awk -v e="\$E" '
      \$0 ~ /^#/ { next }
      NF >= 14 { q=\$3; eval=\$13+0; if (eval <= e) print q }
    ' "${tblout}" | sort -u > "${read_ids}"

    # Filter FASTA to reads with hits
    if [[ -s "${read_ids}" ]]; then
        seqkit grep -f "${read_ids}" "${sampleID}.merged.fa" > "${filtered_fa}"
    else
        : > "${filtered_fa}"
    fi
    """
}
