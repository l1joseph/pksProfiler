nextflow.enable.dsl = 2

// ---------------- Parameters ----------------
params.sample = null

params.input_data_type = "bam"         // bam | fastq
params.pks_taxa = false   // set true to run krakenuniq/bracken on pks-island reads
params.save_intermediates = false // publish extracted/filtered/host-depleted FASTQs

params.profiling_method = "bowtie2" // bowtie2 | hmm | both
params.hmm_evalue       = 1e-10
params.hmm_chunking = false
params.entero_filter = true
params.hmm_model        = "${projectDir}/ref/hmm/clb_all_dna.hmm"

// MAG branch (metagenome mode only)
params.enable_mags       = false
params.gtdbtk_db         = null
params.checkm2_db        = null
params.clb_protein_hmm   = "${projectDir}/ref/hmm/clb_all_protein.hmm"
params.hmm_protein_evalue = 1e-5
params.bracken_read_length = null // Must match a read length supported by the selected Bracken database.

// profile taxa that map to the pks island:
params.pks_shift = 2193827          // island start in E. coli genome coords
params.pks_island_len = 50767       // island length (0..50767 in your file)
params.pks_contig = 'NC_017628.1'   // contig name in BAM

// Output directories
params.outdir = "${launchDir}/results"

params.unmapped_bam_dir = "${params.outdir}/unmapped_reads"
params.mapped_reads_dir = "${params.outdir}/host_depleted_reads"
params.pks_dir = "${params.outdir}/pks_per_sample"
params.pks_summary_dir = "${params.outdir}/pks_summary"
params.pks_counts_dir = "${params.pks_summary_dir}/gene_counts"
params.pks_coverage_plots_dir = "${params.pks_summary_dir}/coverage_plots"
params.pks_taxonomy_dir = "${params.pks_summary_dir}/taxonomy"
params.pks_taxonomy_plots_dir = "${params.pks_taxonomy_dir}/plots"
params.pks_qc_dir = "${params.pks_summary_dir}/qc"
params.pks_mag_dir = "${params.pks_summary_dir}/mag"

// Databases and refs [CHANGE THIS]
params.hg38_db      = null
params.t2t_phix_db  = null
params.pangenome_db = null
params.adapters     = "${projectDir}/ref/known_adapters.fna"
params.kraken_db= null


// PKS + E. coli annotation
params.pks_genome            = "${projectDir}/indices/GCF_000025745.1/GCF_000025745.1_ASM2574v1_genomic"
params.pks_genome_annotation = "${projectDir}/ref/annotations/IHE3034.clbA-clbS.gff"
params.pks_cytoband          = "${projectDir}/indices/GCF_000025745.1/genomic_pks.txt"

// Envs
params.samtools_env = "${projectDir}/conda_envs/samtools_env.yml"
params.fastp_env = "${projectDir}/conda_envs/fastp_env.yml"
params.minimap2_env = "${projectDir}/conda_envs/minimap2_env.yml"
params.pks_align_env = "${projectDir}/conda_envs/pks_align_env.yml"
params.pks_hmm_env = "${projectDir}/conda_envs/pks_hmm_env.yml"
params.krakenuniq_bracken_env = "${projectDir}/conda_envs/krakenUniq_bracken_env.yml"
params.scripts = "${projectDir}/scripts"

// ---------------- Modules ----------------
include { extractReads } from './Modules/extract_reads.nf'
include { filterReads } from './Modules/filter_reads.nf'
include { mapReads } from './Modules/map_reads.nf'
include { pksProfiler_align as pksProfilerAlign } from './Modules/pksProfiler_align.nf'
include { pksProfiler_hmm as pksProfilerHMM } from './Modules/pksProfiler_hmm.nf'
include { plotPKS; masterTableAlign; masterTableHMM; masterQCSummary } from './Modules/plotting.nf'
include { filterEnterobacteriaceae; extractPksIslandReads; Bracken; process_bracken as combinePKSTaxa; combineClbTaxonomySupport } from './Modules/pks_taxa.nf'
include { plotBrackenTaxa as plotPKSTaxa } from './Modules/plot_bracken_taxa.nf'
include { pksMAG } from './Modules/pks_mag.nf'

// ---------------- Workflow ----------------
workflow {

	// ---------- Required input validation ----------
    if (!params.sample) {
        exit 1, "Missing required parameter: --sample"
    }

    if (!params.hg38_db) {
        exit 1, "Missing required parameter: --hg38_db"
    }

    if (!params.t2t_phix_db) {
        exit 1, "Missing required parameter: --t2t_phix_db"
    }

    if (params.pks_taxa && !params.kraken_db) {
        exit 1, "Taxonomic profiling requires: --kraken_db"
    }
	if (params.pks_taxa && !params.bracken_read_length) {
	    exit 1, "Taxonomic profiling requires: --bracken_read_length"
	}
	if (params.pks_taxa && (
        !(params.bracken_read_length.toString() ==~ /^[0-9]+$/) ||
        params.bracken_read_length.toString().toInteger() <= 0
    )) {
	    exit 1, "--bracken_read_length must be a positive integer"
	}

    def hmm_evalue_text = params.hmm_evalue.toString()

    if (!(
        hmm_evalue_text ==~
        /^(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?$/
    )) {
        exit 1, "--hmm_evalue must be a positive number, for example 1e-10"
    }

    if (hmm_evalue_text.toDouble() <= 0) {
        exit 1, "--hmm_evalue must be greater than zero"
    }

    if (params.enable_mags) {
        if (!params.gtdbtk_db) {
            exit 1, "--enable_mags requires --gtdbtk_db"
        }
        if (!params.checkm2_db) {
            exit 1, "--enable_mags requires --checkm2_db"
        }
        if (!file(params.clb_protein_hmm).exists()) {
            exit 1, "--clb_protein_hmm not found: ${params.clb_protein_hmm}\n" +
                    "Build it first: bash scripts/build_clb_protein_hmm.sh"
        }
    }

    // ---------- STEP 1: Inputs + filtering ----------
    if (!(params.input_data_type in ["bam", "fastq"])) {
        exit 1, "Unknown --input_data_type: ${params.input_data_type}. Supported: bam, fastq"
    }

    def required_sample_columns = params.input_data_type == "bam" ?
        ["patient", "bam"] :
        ["patient", "fastq1"]

    def sample_sheet = channel
        .fromPath(params.sample, checkIfExists: true)
        .splitCsv(header: true)
        .collect()
        .flatMap { rows ->
            if (!rows) {
                error "Sample sheet contains no samples: ${params.sample}"
            }

            def observed_columns = rows[0].keySet()
            def missing_columns = required_sample_columns.findAll { column ->
                !(column in observed_columns)
            }

            if (missing_columns) {
                error "Sample sheet is missing required column(s) for ${params.input_data_type} input: ${missing_columns.join(', ')}"
            }

            def observed_ids = new HashSet()
            def duplicate_ids = new TreeSet()

            rows.eachWithIndex { row, index ->
                def row_number = index + 2
                def sample_id = row.patient?.toString()?.trim()

                if (!sample_id) {
                    error "Sample sheet row ${row_number} has an empty patient value"
                }

                if (!(sample_id ==~ /^[A-Za-z0-9][A-Za-z0-9._-]*$/)) {
                    error "Invalid patient value '${sample_id}' on row ${row_number}. Use only letters, numbers, periods, underscores, and hyphens; the first character must be alphanumeric."
                }

                if (!observed_ids.add(sample_id)) {
                    duplicate_ids.add(sample_id)
                }

                required_sample_columns
                    .findAll { column -> column != "patient" }
                    .each { column ->
                        if (!row[column]?.toString()?.trim()) {
                            error "Sample sheet row ${row_number} has an empty ${column} value"
                        }
                    }
            }

            if (duplicate_ids) {
                error "Sample identifiers must be unique. Duplicate patient value(s): ${duplicate_ids.join(', ')}"
            }

            rows
        }

    def QC_FRAGMENTS = channel.empty()

    if (params.input_data_type == "bam") {

		// Expect columns: patient,bam
		sample_sheet = sample_sheet.map { row ->
		    tuple(row.patient, file(row.bam, checkIfExists: true))
		}

        EXTRACT_OUT = extractReads(sample_sheet)

        EXTRACT_OUT.reads
            .map { sampleID, reads -> tuple(sampleID, [reads]) }
            .set { READS_TO_FILTER }

        QC_FRAGMENTS = QC_FRAGMENTS.mix(
            EXTRACT_OUT.qc.map { _sampleID, qc_file -> qc_file }
        )

    } else if (params.input_data_type == "fastq") {

        def sample_sheet_fastq = sample_sheet
            .map { row ->
                def fastq_files = [
                    file(row.fastq1, checkIfExists: true)
                ]

                def fastq2 = row.fastq2?.toString()?.trim()
                if (fastq2) {
                    fastq_files << file(fastq2, checkIfExists: true)
                }

                tuple(
                    row.patient,
                    fastq_files
                )
            }

        sample_sheet_fastq.set { READS_TO_FILTER }

    }

    FILTER_OUT = filterReads(READS_TO_FILTER)

    FILTER_OUT.reads
        .set { FILTERED_UNMAPPED_READS }

    QC_FRAGMENTS = QC_FRAGMENTS.mix(
        FILTER_OUT.qc.map { _sampleID, qc_file -> qc_file }
    )

	// ---------- STEP 1b: Host read depletion ----------
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

	// ---------- STEP 2: Profiling ----------
    def valid_methods = ["bowtie2", "hmm", "both"]

    if (!(params.profiling_method in valid_methods)) {
        exit 1, "Unknown --profiling_method: ${params.profiling_method}. Supported: bowtie2, hmm, both"
    }

    def do_align = params.profiling_method in ["bowtie2", "both"]
    def do_hmm   = params.profiling_method in ["hmm", "both"]

    if (do_align) {
        ALIGN_OUT = pksProfilerAlign(PROFILING_READS)

        ALIGN_OUT.profile
            .set { PKS_ALIGN_OUT }

        QC_FRAGMENTS = QC_FRAGMENTS.mix(
            ALIGN_OUT.qc.map { _sampleID, qc_file -> qc_file }
        )
    }
    if (do_hmm) {
        HMM_OUT = pksProfilerHMM(PROFILING_READS)

        HMM_OUT.profile
            .set { PKS_HMM_OUT }

        QC_FRAGMENTS = QC_FRAGMENTS.mix(
            HMM_OUT.qc.map { _sampleID, qc_file -> qc_file }
        )
    }

    // ---------- STEP 3: Plotting (align only) ----------
    if (do_align) {
        PKS_ALIGN_OUT
            .map { output -> output[2]  }   // bedgraph
            .set { COVERAGE_BEDGRAPH }

        plotPKS(COVERAGE_BEDGRAPH)
    }

	// ---------- STEP 3b: Optional PKS-island taxa profiling (align only) ----------
    if (do_align && params.pks_taxa) {
        PKS_ALIGN_OUT
			.map { sampleID, _covtxt, _bedgraph, _counts, bam, bai, _sam ->
			    tuple(sampleID, bam, bai)
			}
            .set { PKS_BAM_FOR_TAXA }

        extractPksIslandReads(PKS_BAM_FOR_TAXA)
            .set { PKS_ISLAND_FASTQ }

		Bracken(PKS_ISLAND_FASTQ).set { BRACKEN_PER_SAMPLE }

		BRACKEN_PER_SAMPLE
		.map { _sampleID, _report, _classified, _unclassified, _brG, _brS, _gk, _sk, _gmpa, _smpa, speciesSupport -> speciesSupport }
		.collect()
		.set { CLB_SPECIES_SUPPORT_FILES }

		def combine_clb_support_script = file(
		    "${params.scripts}/combine_clb_species_support.py",
		    checkIfExists: true
		)

		combineClbTaxonomySupport(
		    CLB_SPECIES_SUPPORT_FILES,
		    combine_clb_support_script
		)

		BRACKEN_PER_SAMPLE
		.map { sampleID, _kreport, _classified, _unclassified, brG, brS, _gk, _sk, _gmpa, _smpa, _speciesMatrix ->
			    tuple(sampleID, brG, brS)
		}
	    .set { BRACKEN_GS_REPORTS }

		plotPKSTaxa(BRACKEN_GS_REPORTS)

	
		BRACKEN_PER_SAMPLE
		.map { _sampleID, _report, _classified, _unclassified, _brG, _brS, _gk, _sk, gmpa, smpa, _speciesMatrix ->
		    [gmpa, smpa]
		}
	   .flatten()
       .collect()
       .set { BRACKEN_MPA_FILES }

		combinePKSTaxa(BRACKEN_MPA_FILES)
		

    }

    // ---------- STEP 3c: MAG branch (metagenome mode + --enable_mags) ----------
    if (params.enable_mags) {
        pksMAG(MAPPED_READS)
    }

    // ---------- STEP 4: Master tables ----------
    if (do_align) {
        PKS_ALIGN_OUT
            .map { output -> output[3] }    // counts.txt
            .collect()
            .set { ALIGN_COUNT_FILES }

        masterTableAlign(ALIGN_COUNT_FILES)
    }

    if (do_hmm) {
        PKS_HMM_OUT
            .map { output -> output[3] }    // hmm_counts.tsv
            .collect()
            .set { HMM_COUNT_FILES }

        masterTableHMM(HMM_COUNT_FILES)
    }

    // ---------- STEP 5: Cohort QC ----------
    QC_FRAGMENTS
        .collect()
        .set { QC_FRAGMENT_FILES }

    def qc_summary_script = file(
        "${params.scripts}/build_qc_summary.py",
        checkIfExists: true
    )

    masterQCSummary(QC_FRAGMENT_FILES, qc_summary_script)
}
