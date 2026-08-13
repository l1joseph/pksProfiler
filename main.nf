nextflow.enable.dsl = 2

// ---------------- Parameters ----------------
params.sample = "/tscc/nfs/home/amabbasi/restricted/microbiome_pipeline/sample.csv"

params.input_context   = "host"  // metagenome | host
params.input_data_type = "bam"         // bam | fastq
params.pks_taxa = true   // set true to run krakenuniq/bracken on pks-island reads

params.profiling_method = "both"    // bowtie2 | hmm | both
params.hmm_evalue       = 1e-10
params.hmm_model        = "${projectDir}/ref/hmm/clb_all_dna.hmm"

// profile taxa that map to the pks island:
params.pks_shift = 2193827          // island start in E. coli genome coords
params.pks_island_len = 50767       // island length (0..50767 in your file)
params.pks_contig = 'NC_017628.1'   // contig name in BAM

// Output directories
params.unmapped_bam_dir = "${projectDir}/RESULTS/UNMAPPED_BAM"
params.mapped_reads_dir = "${projectDir}/RESULTS/MAPPED_READS"
params.pks_dir = "${projectDir}/RESULTS/PKS_PER_SAMPLE"
params.pks_summary_dir = "${projectDir}/RESULTS/PKS_SUMMARY"

// Databases and refs [CHANGE THIS]
params.hg38_db      = "/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/human-GRC-db.mmi"
params.t2t_phix_db  = "/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/human-GCA-phix-db.mmi"
params.pangenome_db = "/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/pangenome_mmi"
params.adapters     = "${projectDir}/ref/known_adapters.fna"
params.kraken_db="/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/dbs/krakenUniq_8_8_2023"


// PKS + E. coli annotation
params.pks_genome            = "${projectDir}/indices/GCF_000025745.1/GCF_000025745.1_ASM2574v1_genomic"
params.pks_genome_annotation = "${projectDir}/indices/GCF_000025745.1/genomic.gff"
params.pks_cytoband          = "${projectDir}/indices/GCF_000025745.1/genomic_pks.txt"
params.ecoli_cytoband        = "${projectDir}/indices/GCF_000025745.1/genomic_ecoli.txt"

// Envs
params.samtools_env  = "./conda_envs/samtools_env.yml"
params.fastp_env     = "./conda_envs/fastp_env.yml"
params.minimap2_env  = "./conda_envs/minimap2_env.yml"
params.pks_align_env = "./conda_envs/pks_align_env.yml"
params.pks_hmm_env   = "./conda_envs/pks_hmm_env.yml"
params.krakenuniq_bracken_env = "./conda_envs/krakenUniq_bracken_env.yml"
params.krakentools_pack ="/tscc/projects/ps-lalexandrov/shared/CMPipeline_nextflow/packages/KrakenTools"

params.scripts = "${projectDir}/scripts"

// ---------------- Modules ----------------
include { extractReads } from './Modules/extract_reads.nf'
include { filterReads } from './Modules/filter_reads.nf'
include { mapReads } from './Modules/map_reads.nf'
include { mapReads as mapReads_metagenome } from './Modules/map_reads_metagenome.nf'
include { pksProfiler_align as pksProfilerAlign } from './Modules/pksProfiler_align.nf'
include { pksProfiler_hmm as pksProfilerHMM } from './Modules/pksProfiler_hmm.nf'
include { plotPKS; plotGenome; masterTableAlign; masterTableHMM } from './Modules/plotting.nf'
include { extractPksIslandReads; Bracken; process_bracken as combinePKSTaxa } from './Modules/pks_taxa.nf'
include { plotBrackenTaxa as plotPKSTaxa } from './Modules/plot_bracken_taxa.nf'

// ---------------- Workflow ----------------
workflow {

    // ---------- STEP 1: Inputs + filtering ----------
    def sample_sheet = Channel
        .fromPath(params.sample)
        .splitCsv(header: true)

    if (params.input_data_type == "bam") {

        // Expect columns: patient,bam
        sample_sheet = sample_sheet.map { row -> row.subMap('patient', 'bam') }

        extractReads(sample_sheet).set { UNMAPPED_READS }

        UNMAPPED_READS
            .multiMap { sampleID, r1, r2 ->
                whole:     tuple(sampleID, r1, r2)
                path_only: tuple(r1, r2)
            }
            .set { UNMAPPED_READS_MULTI }

        filterReads(UNMAPPED_READS_MULTI.whole)
            .set { FILTERED_UNMAPPED_READS }

    } else if (params.input_data_type == "fastq") {

        def sample_sheet_fastq = sample_sheet
            .map { row -> row.subMap('patient', 'fastq1', 'fastq2') }
            .map { row -> tuple(row.patient, file(row.fastq1), file(row.fastq2)) }

        filterReads(sample_sheet_fastq)
            .set { FILTERED_UNMAPPED_READS }

    } else {
        exit 1, "Unknown --input_data_type: ${params.input_data_type}. Supported: bam, fastq"
    }

    FILTERED_UNMAPPED_READS
        .multiMap { sampleID, r1, r2 ->
            whole:     tuple(sampleID, r1, r2)
            path_only: tuple(r1, r2)
        }
        .set { FILTERED_UNMAPPED_READS_MULTI }

    // ---------- STEP 1b: Mapping ----------
   // gather the list of pangenome .mmi files
    def mmiFiles = []
    def dir = new File("${params.pangenome_db}")
    dir.eachFileRecurse (groovy.io.FileType.FILES) { file ->
        if (file.name.endsWith('.mmi')) {
            mmiFiles << file
        }
    }

    if (params.input_context == "metagenome") {
        mapReads_metagenome(FILTERED_UNMAPPED_READS_MULTI.whole)
            .set { MAPPED_READS }

    } else if (params.input_context == "host") {
        mapReads(FILTERED_UNMAPPED_READS_MULTI.whole, mmiFiles)
            .set { MAPPED_READS }

    } else {
        exit 1, "Unknown --input_context: ${params.input_context}. Supported: metagenome, host"
    }

    MAPPED_READS
        .multiMap { sampleID, r1Pan, r2Pan ->
            PAN: tuple(sampleID, r1Pan, r2Pan)
        }
        .set { MAPPED_READS_MULTI }

    // ---------- STEP 2: Profiling ----------
    def valid_methods = ["bowtie2", "hmm", "both"]
    if (!(params.profiling_method in valid_methods)) {
        exit 1, "Unknown --profiling_method: ${params.profiling_method}. Supported: bowtie2, hmm, both"
    }

    def do_align = params.profiling_method in ["bowtie2", "both"]
    def do_hmm   = params.profiling_method in ["hmm", "both"]

    if (do_align) {
        pksProfilerAlign(MAPPED_READS_MULTI.PAN)
            .set { PKS_ALIGN_OUT }
    }
    if (do_hmm) {
        pksProfilerHMM(MAPPED_READS_MULTI.PAN)
            .set { PKS_HMM_OUT }
    }

    // ---------- STEP 3: Plotting (align only) ----------
    if (do_align) {
        PKS_ALIGN_OUT
            .map { it[2] }   // bedgraph
            .set { COVERAGE_BEDGRAPH }

        plotPKS(COVERAGE_BEDGRAPH)
        plotGenome(COVERAGE_BEDGRAPH)
    }

	// ---------- STEP 3b: Optional PKS-island taxa profiling (align only) ----------
    if (do_align && params.pks_taxa) {
        PKS_ALIGN_OUT
            .map { sampleID, covtxt, bedgraph, counts, bam, bai, sam ->
                tuple(sampleID, bam, bai)
            }
            .set { PKS_BAM_FOR_TAXA }

        extractPksIslandReads(PKS_BAM_FOR_TAXA)
            .set { PKS_ISLAND_FASTQ }

		Bracken(PKS_ISLAND_FASTQ).set { BRACKEN_PER_SAMPLE }

		BRACKEN_PER_SAMPLE
		  .map { sampleID, kreport, classified, unclassified, brG, brS, gk, sk, gmpa, smpa ->
	      tuple(sampleID, brG, brS)
		}
	    .set { BRACKEN_GS_REPORTS }

		plotPKSTaxa(BRACKEN_GS_REPORTS)

	
		BRACKEN_PER_SAMPLE
		  .map { sampleID, report, classified, unclassified, brG, brS, gk, sk, gmpa, smpa ->
	      [ gmpa, smpa ]
		}
	   .flatten()
       .collect()
       .set { BRACKEN_MPA_FILES }

		combinePKSTaxa(BRACKEN_MPA_FILES)
		

    }

    // ---------- STEP 4: Master tables ----------
    if (do_align) {
        PKS_ALIGN_OUT
            .map { it[3] }   // counts.txt
            .collect()
            .set { ALIGN_COUNT_FILES }

        masterTableAlign(ALIGN_COUNT_FILES)
    }

    if (do_hmm) {
        PKS_HMM_OUT
            .map { it[3] }   // hmm_counts.tsv
            .collect()
            .set { HMM_COUNT_FILES }

        masterTableHMM(HMM_COUNT_FILES)
    }
}





