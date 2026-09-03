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
    tuple val(sampleID), path("bins/bin.*.fa"), emit: bins

    script:
    """
    set -euo pipefail
    mkdir -p bins
    metabat2 -i ${contigs} -a ${depth} -o bins/bin -t ${task.cpus} --unbinned || true
    # Ensure at least one bin exists; if metabat2 produced nothing, treat all
    # contigs as a single "bin.0.fa" so downstream processes have input.
    if ! ls bins/bin.*.fa 2>/dev/null | grep -q .; then
        cp ${contigs} bins/bin.0.fa
    fi
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
    for b in ${bins}; do cp \$b bin_input/; done
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
    for b in ${bins}; do cp \$b bin_input/; done
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

// ─── Prodigal gene prediction (per bin) ───────────────────────────────────────

process prodigalPredict {
    label 'mag_hmm'
    scratch true
    publishDir { "${params.outdir}/pks_summary/mag/${sampleID}" }, mode: 'copy',
        saveAs: { "${binID}.faa" }
    conda "${projectDir}/conda_envs/prodigal_env.yml"

    input:
    tuple val(sampleID), val(binID), path(bin_fa)

    output:
    tuple val(sampleID), val(binID), path("${binID}.faa"), emit: proteins

    script:
    """
    set -euo pipefail
    prodigal -i ${bin_fa} -a ${binID}.faa -p meta -f gff -q
    """
}

// ─── hmmsearch vs colibactin protein HMM (per bin) ────────────────────────────

process hmmsearchClb {
    label 'mag_hmm'
    scratch true
    conda "${projectDir}/conda_envs/pks_hmm_env.yml"

    input:
    tuple val(sampleID), val(binID), path(proteins)

    output:
    tuple val(sampleID), val(binID), path("${binID}.tblout"), emit: tblout

    script:
    """
    set -euo pipefail
    hmmsearch --cpu ${task.cpus} \
        --tblout ${binID}.tblout \
        -E ${params.hmm_protein_evalue} \
        ${params.clb_protein_hmm} \
        ${proteins}
    """
}

// ─── Prokka annotation (pks+ bins only) ───────────────────────────────────────

process prokkaAnnotate {
    label 'mag_hmm'
    scratch true
    publishDir { "${params.outdir}/pks_summary/mag/${sampleID}" }, mode: 'copy',
        saveAs: { "${binID}.gff" }
    conda "${projectDir}/conda_envs/prokka_env.yml"

    input:
    tuple val(sampleID), val(binID), path(bin_fa)

    output:
    tuple val(sampleID), val(binID), path("prokka_out/${binID}.gff"), emit: gff

    script:
    """
    set -euo pipefail
    prokka --outdir prokka_out --prefix ${binID} \
        --metagenome --cpus ${task.cpus} --force --quiet \
        ${bin_fa}
    """
}

// ─── Genomic context extraction (pks+ bins only) ──────────────────────────────

process extractGenomicContext {
    label 'mag_hmm'
    scratch true
    publishDir { "${params.outdir}/pks_summary/mag/${sampleID}" }, mode: 'copy',
        saveAs: { "${binID}.context.tsv" }
    conda "${projectDir}/conda_envs/pks_hmm_env.yml"

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
    publishDir "${params.outdir}/pks_summary/mag", mode: 'copy'
    conda "${projectDir}/conda_envs/pks_hmm_env.yml"

    input:
    tuple val(sampleID), path(checkm2_report), path(gtdbtk_summary),
          path(tblouts), path(contexts)

    output:
    path "${sampleID}.pks_mag_summary.tsv", emit: summary

    script:
    """
    set -euo pipefail
    mkdir -p tblout_dir context_dir
    for f in ${tblouts}; do cp \$f tblout_dir/; done
    for f in ${contexts}; do cp \$f context_dir/; done
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

    // 7. Prodigal gene prediction (per bin)
    prodigalPredict(bins_flat_ch)

    // 8. hmmsearch vs clb protein HMM (per bin)
    hmmsearchClb(prodigalPredict.out.proteins)

    // 9. Filter to pks+ bins (≥1 passing hit in tblout)
    pks_pos_tblout_ch = hmmsearchClb.out.tblout
        .filter { sampleID, binID, tblout ->
            tblout.readLines().any { line -> !line.startsWith('#') && !line.trim().isEmpty() }
        }

    // 10. Prokka annotation (pks+ bins only)
    pks_pos_fa_ch = pks_pos_tblout_ch
        .map { sampleID, binID, tblout -> tuple(sampleID, binID) }
        .join(bins_flat_ch.map { sampleID, binID, fa -> tuple(sampleID, binID, fa) }, by: [0, 1])
    prokkaAnnotate(pks_pos_fa_ch)

    // 11. Genomic context (gff + tblout joined per bin)
    gff_tblout_ch = prokkaAnnotate.out.gff
        .join(pks_pos_tblout_ch, by: [0, 1])
    extractGenomicContext(gff_tblout_ch)

    // 12. Aggregate per sample and build summary
    tblouts_per_sample_ch = hmmsearchClb.out.tblout
        .map { sampleID, binID, tblout -> tuple(sampleID, tblout) }
        .groupTuple(by: 0)

    contexts_per_sample_ch = extractGenomicContext.out.context
        .map { sampleID, binID, ctx -> tuple(sampleID, ctx) }
        .groupTuple(by: 0)

    summary_input_ch = checkm2Predict.out.report
        .join(gtdbtkClassify.out.summary, by: 0)
        .join(tblouts_per_sample_ch, by: 0)
        .join(contexts_per_sample_ch, by: 0)

    magSummaryTable(summary_input_ch)
}
