process extractPksIslandReads {
  label 'process_medium'
  scratch true
  publishDir "${params.pks_dir}", mode: 'copy'
  conda "${params.pks_align_env}"

  input:
  tuple val(sampleID), path(bam), path(bai)

  output:
  tuple val(sampleID),
        path("${sampleID}.pks.fastq.gz"),
        path("${sampleID}.read_clb_gene.tsv")

  script:
  """
  set -euo pipefail

  CHR="${params.pks_contig}"
  START=\$(( ${params.pks_shift} - 1 ))
  END=\$(( ${params.pks_shift} + ${params.pks_island_len} ))

  printf "%s\\t%s\\t%s\\n" \
    "\$CHR" "\$START" "\$END" \
    > pks_island.bed

  # Extract alignments overlapping the complete pks island
  samtools view \
    -b \
    -L pks_island.bed \
    "${bam}" \
    > "${sampleID}.pks.bam"

  # Convert clb gene coordinates from GFF to BED
  awk -F '\\t' '
    BEGIN {
      OFS="\\t"
    }

    \$0 !~ /^#/ && \$3 == "gene" {
      gene=""
      n=split(\$9, attributes, ";")

      for (i=1; i<=n; i++) {
        if (attributes[i] ~ /^Name=/) {
          sub(/^Name=/, "", attributes[i])
          gene=attributes[i]
        }
      }

      if (gene != "") {
        print \$1, \$4-1, \$5, gene
      }
    }
  ' "${params.pks_genome_annotation}" > clb_genes.bed

  # Assign each read/template to the clb gene with the greatest
  # total aligned-base overlap
  bedtools bamtobed \
    -i "${sampleID}.pks.bam" |
    bedtools intersect \
      -a - \
      -b clb_genes.bed \
      -wo |
    awk '
      BEGIN {
        OFS="\\t"
      }

      {
        read_id=\$4
        gene=\$10
        overlap=\$11+0
        total[read_id SUBSEP gene] += overlap
      }

      END {
        for (key in total) {
          split(key, fields, SUBSEP)

          read_id=fields[1]
          gene=fields[2]
          overlap=total[key]

          if (!(read_id in best_overlap) ||
              overlap > best_overlap[read_id] ||
              (overlap == best_overlap[read_id] &&
               gene < best_gene[read_id])) {
            best_overlap[read_id]=overlap
            best_gene[read_id]=gene
          }
        }

        for (read_id in best_gene) {
          print read_id, best_gene[read_id], best_overlap[read_id]
        }
      }
    ' > "${sampleID}.read_clb_gene.tmp.tsv"

  printf "read_id\\tGene\\toverlap_bp\\n" \
    > "${sampleID}.read_clb_gene.tsv"

  sort -k1,1 "${sampleID}.read_clb_gene.tmp.tsv" \
    >> "${sampleID}.read_clb_gene.tsv"

  rm -f "${sampleID}.read_clb_gene.tmp.tsv"

  # Convert overlapping alignments to a single FASTQ stream
  samtools fastq "${sampleID}.pks.bam" |
    gzip -c > "${sampleID}.pks.fastq.gz"
  """
}


process Bracken {
  scratch true
  label 'process_high_disk'
  publishDir "${params.pks_dir}", mode: 'copy'
  conda "${params.krakenuniq_bracken_env}"

  input:
  tuple val(sampleID), path(fastq_gz), path(read_gene_tsv)

  output:
  tuple val(sampleID),
        path("${sampleID}.krakenuniq.report.txt"),
        path("${sampleID}.classified.fasta"),
        path("${sampleID}.unclassified.fasta"),
        path("${sampleID}.bracken.G.report.txt"),
        path("${sampleID}.bracken.S.report.txt"),
        path("${sampleID}.bracken.G.krakenreport.txt"),
        path("${sampleID}.bracken.S.krakenreport.txt"),
        path("${sampleID}.bracken.G.mpa.krakenreport.txt"),
        path("${sampleID}.bracken.S.mpa.krakenreport.txt"),
        path("${sampleID}.clb_species_support.tsv")

  script:
  """
  set -euo pipefail

  REPORT="${sampleID}.krakenuniq.report.txt"
  OUTPUT="${sampleID}.krakenuniq.output.txt"
  CLASSIFIED="${sampleID}.classified.fasta"
  UNCLASSIFIED="${sampleID}.unclassified.fasta"
  SPECIES_MATRIX="${sampleID}.clb_species_support.tsv"

  # Decompress PKS reads for KrakenUniq
  zcat "${fastq_gz}" > "${sampleID}.pks.fastq"

  # A valid sample can contain no PKS reads
  if [[ ! -s "${sampleID}.pks.fastq" ]]; then
    echo "No PKS reads for ${sampleID}; writing empty outputs."

    : > "\$REPORT"
    : > "\$OUTPUT"
    : > "\$CLASSIFIED"
    : > "\$UNCLASSIFIED"

    for lvl in G S; do
      : > "${sampleID}.bracken.\${lvl}.report.txt"
      : > "${sampleID}.bracken.\${lvl}.krakenreport.txt"
      : > "${sampleID}.bracken.\${lvl}.mpa.krakenreport.txt"
    done

    # Write a valid header-only species-by-clb matrix
    {
      printf "Species\\tTaxID"

      for gene in {A..S}; do
        printf "\\tclb%s" "\$gene"
      done

      printf "\\tTotal\\n"
    } > "\$SPECIES_MATRIX"

    exit 0
  fi

  krakenuniq \
    --db "${params.kraken_db}" \
    --threads "${task.cpus}" \
    --report-file "\$REPORT" \
    --output "\$OUTPUT" \
    --classified-out "\$CLASSIFIED" \
    --unclassified-out "\$UNCLASSIFIED" \
    "${sampleID}.pks.fastq"

  # Direct per-read KrakenUniq assignments retain the connection
  # between taxon and clb gene. Bracken outputs remain separate because
  # Bracken estimates aggregate abundance and has no per-read identity.
  python "${params.scripts}/build_clb_species_matrix.py" \
    --read-gene "${read_gene_tsv}" \
    --kraken-output "\$OUTPUT" \
    --kraken-db "${params.kraken_db}" \
    --output "\$SPECIES_MATRIX"

  # Apply the same conservative support threshold used by Bracken.
  # Reads classified only above the requested rank do not satisfy this
  # requirement until at least two reads support a target-rank node.
  GENUS_READS=\$(awk -F '\\t' '
    \$8 == "genus" && \$2 ~ /^[0-9]+\$/ {
      sum += \$2
    }

    END {
      print sum+0
    }
  ' "\$REPORT")

  SPECIES_READS=\$(awk -F '\\t' '
    \$8 == "species" && \$2 ~ /^[0-9]+\$/ {
      sum += \$2
    }

    END {
      print sum+0
    }
  ' "\$REPORT")

  for lvl in G S; do
    bracken_output="${sampleID}.bracken.\${lvl}.report.txt"
    bracken_kraken_report="${sampleID}.bracken.\${lvl}.krakenreport.txt"
    bracken_kraken_mpa_report="${sampleID}.bracken.\${lvl}.mpa.krakenreport.txt"

    if [[ "\$lvl" == "G" ]]; then
      LVL_READS="\$GENUS_READS"
    else
      LVL_READS="\$SPECIES_READS"
    fi

    if [[ "\$LVL_READS" -lt 2 ]]; then
      echo "Skipping Bracken level \$lvl: exact-rank reads=\$LVL_READS; threshold=2"

      : > "\$bracken_output"
      : > "\$bracken_kraken_report"
      : > "\$bracken_kraken_mpa_report"

      continue
    fi

    bracken \
      -d "${params.kraken_db}" \
      -i "\$REPORT" \
      -o "\$bracken_output" \
      -w "\$bracken_kraken_report" \
      -r ${params.bracken_read_length} \
      -l "\$lvl" \
      -t 2

    kreport2mpa.py \
      -r "\$bracken_kraken_report" \
      -o "\$bracken_kraken_mpa_report" \
      --display-header
  done
  """
}


process process_bracken {
  scratch true
  publishDir "${params.pks_taxonomy_dir}", mode: 'copy'
  conda "${params.krakenuniq_bracken_env}"

  input:
  path bracken_files

  output:
  tuple path("bracken.genus.mpa.report.txt"),
        path("bracken.species.mpa.report.txt")

  script:
  def genus_files = bracken_files.findAll { file ->
    file.name.endsWith('.G.mpa.krakenreport.txt')
  }

  def species_files = bracken_files.findAll { file ->
    file.name.endsWith('.S.mpa.krakenreport.txt')
  }

  def genus_str = genus_files
    .collect { file -> "\"${file}\"" }
    .join(' ')

  def species_str = species_files
    .collect { file -> "\"${file}\"" }
    .join(' ')

  """
  set -euo pipefail

  if [[ -n "${genus_str}" ]]; then
    combine_mpa.py \
      --input ${genus_str} \
      --output bracken.genus.mpa.report.txt
  else
    echo "No genus files found." \
      > bracken.genus.mpa.report.txt
  fi

  if [[ -n "${species_str}" ]]; then
    combine_mpa.py \
      --input ${species_str} \
      --output bracken.species.mpa.report.txt
  else
    echo "No species files found." \
      > bracken.species.mpa.report.txt
  fi
  """
}


process combineClbTaxonomySupport {
  scratch true
  publishDir "${params.pks_taxonomy_dir}", mode: 'copy'
  conda "${params.krakenuniq_bracken_env}"

  input:
  path species_support_files
  path combine_script

  output:
  path("pks.clb_species_support.tsv")

  script:
  def species_inputs = species_support_files
    .collect { file -> "\"${file}\"" }
    .join(' ')

  """
  set -euo pipefail

  python "${combine_script}" \
    --species-files ${species_inputs} \
    --species-output pks.clb_species_support.tsv
  """
}

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

    extract_kraken_reads.py \\
        -k "\$OUTPUT" \\
        -r "\$REPORT" \\
        -s "${sampleID}.all_reads.fastq" \\
        -t 543 \\
        --include-children \\
        -o "\$FILTERED"

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
