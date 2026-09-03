nextflow.enable.dsl=2

// ─── Assembly ─────────────────────────────────────────────────────────────────

process megahitAssemble {
    label 'mag_assembly'
    scratch true
    publishDir { "${params.outdir}/pks_summary/mag/${sampleID}/assembly" },
        mode: 'copy', enabled: params.save_intermediates
    conda "${projectDir}/conda_envs/megahit_env.yml"

    input:
    tuple val(sampleID), path(reads)

    output:
    tuple val(sampleID), path("${sampleID}.contigs.fa"), emit: contigs

    script:
    def reads_list = reads instanceof List ? reads : [reads]
    def r_flag = reads_list.size() == 2
        ? "-1 ${reads_list[0]} -2 ${reads_list[1]}"
        : "-r ${reads_list[0]}"
    """
    set -euo pipefail
    megahit ${r_flag} -t ${task.cpus} -o megahit_out
    cp megahit_out/final.contigs.fa ${sampleID}.contigs.fa
    """
}

// ─── Read-to-contig alignment ─────────────────────────────────────────────────

process alignToContigs {
    label 'mag_binning'
    scratch true
    conda "${projectDir}/conda_envs/minimap2_env.yml"

    input:
    tuple val(sampleID), path(reads), path(contigs)

    output:
    tuple val(sampleID), path("${sampleID}.contigs.sorted.bam"),
          path("${sampleID}.contigs.sorted.bam.bai"), emit: bam

    script:
    def reads_list = reads instanceof List ? reads : [reads]
    def read_args = reads_list.join(' ')
    """
    set -euo pipefail
    minimap2 -ax sr -t ${task.cpus} ${contigs} ${read_args} \
        | samtools view -F 4 -bS \
        | samtools sort -@ ${task.cpus} -o ${sampleID}.contigs.sorted.bam
    samtools index ${sampleID}.contigs.sorted.bam
    """
}

// ─── Contig coverage depth ────────────────────────────────────────────────────

process jgiContigDepths {
    label 'mag_binning'
    scratch true
    conda "${projectDir}/conda_envs/metabat2_env.yml"

    input:
    tuple val(sampleID), path(bam), path(bai)

    output:
    tuple val(sampleID), path("${sampleID}.depth.txt"), emit: depth

    script:
    """
    set -euo pipefail
    jgi_summarize_bam_contig_depths --outputDepth ${sampleID}.depth.txt ${bam}
    """
}

// ─── MetaBAT2 binning ─────────────────────────────────────────────────────────

process metabat2Bin {
    label 'mag_binning'
    scratch true
    publishDir { "${params.outdir}/pks_summary/mag/${sampleID}/bins" },
        mode: 'copy', enabled: params.save_intermediates
    conda "${projectDir}/conda_envs/metabat2_env.yml"

    input:
    tuple val(sampleID), path(contigs), path(depth)

    output:
    // optional: true — samples with too few contigs for binning emit nothing and
    // are silently dropped from the MAG branch rather than passed through as a
    // fake pseudo-bin.
    tuple val(sampleID), path("bins/bin.*.fa"), emit: bins, optional: true

    script:
    """
    set -euo pipefail
    mkdir -p bins
    metabat2 -i ${contigs} -a ${depth} -o bins/bin -t ${task.cpus} || {
        rc=\$?
        if [[ \${rc} -eq 1 ]]; then
            echo "MetaBAT2 exited 1 (no bins produced; expected for low-coverage assemblies)" >&2
        else
            echo "MetaBAT2 failed with exit code \${rc}" >&2
            exit \${rc}
        fi
    }
    """
}

// ─── CheckM2 quality assessment ───────────────────────────────────────────────

process checkm2Predict {
    label 'mag_binning'
    scratch true
    publishDir { "${params.outdir}/pks_summary/mag/${sampleID}/checkm2" }, mode: 'copy'
    conda "${projectDir}/conda_envs/checkm2_env.yml"

    input:
    tuple val(sampleID), path(bins)

    output:
    tuple val(sampleID), path("checkm2_out/quality_report.tsv"), emit: report

    script:
    """
    set -euo pipefail
    mkdir -p bin_input
    for b in ${bins}; do ln -s \$(realpath \$b) bin_input/; done
    checkm2 predict --threads ${task.cpus} \
        --input bin_input \
        --output-directory checkm2_out \
        --database_path ${params.checkm2_db}
    """
}

// ─── GTDB-Tk taxonomy ─────────────────────────────────────────────────────────

process gtdbtkClassify {
    label 'mag_gtdbtk'
    scratch true
    publishDir { "${params.outdir}/pks_summary/mag/${sampleID}/gtdbtk" }, mode: 'copy'
    conda "${projectDir}/conda_envs/gtdbtk_env.yml"

    input:
    tuple val(sampleID), path(bins)

    output:
    tuple val(sampleID), path("gtdbtk_out/gtdbtk.bac120.summary.tsv"), emit: summary

    script:
    """
    set -euo pipefail
    mkdir -p bin_input gtdbtk_out
    for b in ${bins}; do ln -s \$(realpath \$b) bin_input/; done
    GTDBTK_DATA_PATH=${params.gtdbtk_db} gtdbtk classify_wf \
        --genome_dir bin_input \
        --out_dir gtdbtk_out \
        --cpus ${task.cpus} \
        --extension fa \
        --skip_ani_screen
    # Stub if no bacterial bins (only archaea)
    [[ -f gtdbtk_out/gtdbtk.bac120.summary.tsv ]] || \
        printf "user_genome\tclassification\n" > gtdbtk_out/gtdbtk.bac120.summary.tsv
    """
}

// ─── Prokka annotation (ALL bins — must run before hmmsearch so locus_tags match) ──

process prokkaAnnotate {
    label 'mag_hmm'
    scratch true
    publishDir { "${params.outdir}/pks_summary/mag/${sampleID}" }, mode: 'copy',
        saveAs: { fn -> fn.tokenize('/').last() }
    conda "${projectDir}/conda_envs/prokka_env.yml"

    input:
    tuple val(sampleID), val(binID), path(bin_fa)

    output:
    tuple val(sampleID), val(binID), path("prokka_out/${binID}.gff"), emit: gff
    tuple val(sampleID), val(binID), path("prokka_out/${binID}.faa"), emit: faa_for_hmm

    script:
    """
    set -euo pipefail
    prokka --outdir prokka_out --prefix ${binID} \
        --metagenome --cpus ${task.cpus} --force --quiet \
        ${bin_fa}
    """
}

// ─── hmmsearch vs colibactin protein HMM (per bin, on Prokka proteins) ──────────

process hmmsearchClb {
    label 'mag_hmm'
    scratch true
    conda "${params.pks_hmm_env}"

    input:
    tuple val(sampleID), val(binID), path(proteins)

    output:
    tuple val(sampleID), val(binID), path("${binID}.tblout"), emit: tblout
    tuple val(sampleID), val(binID), path("${binID}.hit_count.txt"), emit: hit_count

    script:
    """
    set -euo pipefail
    hmmsearch --cpu ${task.cpus} \
        --tblout ${binID}.tblout \
        -E ${params.hmm_protein_evalue} \
        ${params.clb_protein_hmm} \
        ${proteins}
    grep -vc '^#' ${binID}.tblout > ${binID}.hit_count.txt 2>/dev/null || echo 0 > ${binID}.hit_count.txt
    """
}

// ─── Genomic context extraction (pks+ bins only) ──────────────────────────────

process extractGenomicContext {
    label 'mag_hmm'
    scratch true
    publishDir { "${params.outdir}/pks_summary/mag/${sampleID}" }, mode: 'copy',
        saveAs: { "${binID}.context.tsv" }
    conda "${params.pks_hmm_env}"

    input:
    tuple val(sampleID), val(binID), path(gff), path(tblout)

    output:
    tuple val(sampleID), val(binID), path("${binID}.context.tsv"), emit: context

    script:
    """
    set -euo pipefail
    python ${projectDir}/scripts/extract_genomic_context.py \
        --gff ${gff} \
        --tblout ${tblout} \
        --evalue ${params.hmm_protein_evalue} \
        --out ${binID}.context.tsv
    """
}

// ─── Per-sample MAG summary table ─────────────────────────────────────────────

process magSummaryTable {
    label 'mag_hmm'
    publishDir "${params.pks_mag_dir}", mode: 'copy'
    conda "${params.pks_hmm_env}"

    input:
    tuple val(sampleID), path(checkm2_report), path(gtdbtk_summary),
          path(tblouts), path(contexts)

    output:
    path "${sampleID}.pks_mag_summary.tsv", emit: summary

    script:
    """
    set -euo pipefail
    mkdir -p tblout_dir context_dir
    for f in ${tblouts}; do ln -s \$(realpath \$f) tblout_dir/; done
    for f in ${contexts}; do ln -s \$(realpath \$f) context_dir/; done
    python ${projectDir}/scripts/build_mag_summary.py \
        --checkm2 ${checkm2_report} \
        --gtdbtk ${gtdbtk_summary} \
        --tblout_dir tblout_dir \
        --context_dir context_dir \
        --sample ${sampleID} \
        --evalue ${params.hmm_protein_evalue} \
        --out ${sampleID}.pks_mag_summary.tsv
    """
}

// ─── pksMAG workflow ──────────────────────────────────────────────────────────

workflow pksMAG {
    take:
    reads  // tuple val(sampleID), path(reads)

    main:

    // 1. Assembly
    megahitAssemble(reads)

    // 2. Align reads to contigs for depth estimation
    reads_contigs_ch = reads.join(megahitAssemble.out.contigs, by: 0)
    alignToContigs(reads_contigs_ch)

    // 3. Contig coverage depth
    jgiContigDepths(alignToContigs.out.bam)

    // 4. Bin contigs
    contigs_depth_ch = megahitAssemble.out.contigs
        .join(jgiContigDepths.out.depth, by: 0)
    metabat2Bin(contigs_depth_ch)

    // Fan-out: one entry per bin → tuple(sampleID, binID, bin_fa)
    bins_flat_ch = metabat2Bin.out.bins
        .transpose()
        .map { sampleID, bin_fa ->
            def binID = bin_fa.baseName
            tuple(sampleID, binID, bin_fa)
        }

    // 5. CheckM2 quality (all bins together per sample)
    checkm2Predict(metabat2Bin.out.bins)

    // 6. GTDB-Tk taxonomy (all bins together per sample)
    gtdbtkClassify(metabat2Bin.out.bins)

    // 7. Prokka annotation (ALL bins — before hmmsearch so locus_tags are consistent)
    prokkaAnnotate(bins_flat_ch)

    // 8. hmmsearch vs clb protein HMM (per bin, on prokka proteins)
    hmmsearchClb(prokkaAnnotate.out.faa_for_hmm)

    // 9. Filter to pks+ bins using hit_count (avoids reading tblouts in the driver JVM)
    pks_pos_tblout_ch = hmmsearchClb.out.tblout
        .join(hmmsearchClb.out.hit_count, by: [0, 1])
        .filter { sampleID, binID, tblout, hit_count ->
            hit_count.text.trim() as Integer > 0
        }
        .map { sampleID, binID, tblout, hit_count -> tuple(sampleID, binID, tblout) }

    // 10. Genomic context (prokka GFF + tblout joined per pks+ bin; locus_tags now match)
    extractGenomicContext(prokkaAnnotate.out.gff.join(pks_pos_tblout_ch, by: [0, 1]))

    // 11. Aggregate per sample and build summary
    tblouts_per_sample_ch = hmmsearchClb.out.tblout
        .map { sampleID, binID, tblout -> tuple(sampleID, tblout) }
        .groupTuple(by: 0)

    contexts_per_sample_ch = extractGenomicContext.out.context
        .map { sampleID, binID, ctx -> tuple(sampleID, ctx) }
        .groupTuple(by: 0)

    // remainder: true keeps all samples even when contexts_per_sample_ch has no
    // entry (i.e. zero pks+ bins); null → [] so the script receives an empty file list.
    summary_input_ch = checkm2Predict.out.report
        .join(gtdbtkClassify.out.summary, by: 0)
        .join(tblouts_per_sample_ch, by: 0)
        .join(contexts_per_sample_ch, by: 0, remainder: true)
        .map { sampleID, checkm2, gtdbtk, tblouts, contexts ->
            tuple(sampleID, checkm2, gtdbtk, tblouts, contexts ?: [])
        }

    magSummaryTable(summary_input_ch)
}
