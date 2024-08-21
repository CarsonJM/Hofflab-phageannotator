/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// FUNCTIONS: Local functions
//
include { getWorkDirs; rmEmptyFastAs; rmEmptyFastQs } from '../../lib/functions.nf'

//
// SUBWORKFLOWS: Consisting of a mix of local and nf-core/modules
//
include { FASTQ_HOSTREMOVAL_BOWTIE2         } from '../../subworkflows/local/fastq_hostremoval_bowtie2/main'
include { FASTQ_VIRUSENRICHMENT_VIROMEQC    } from '../../subworkflows/local/fastq_virusenrichment_viromeqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// MODULES: Installed directly from nf-core/modules
//
include { CAT_FASTQ as CAT_RUNMERGE         } from '../../modules/nf-core/cat/fastq/main'
include { FASTP                             } from '../../modules/nf-core/fastp/main'
include { FASTQC as FASTQC_RAW              } from '../../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_PREPROCESSED     } from '../../modules/nf-core/fastqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN READ PREPROCESSING WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow FASTQ_READPREPROCESSING_FASTQ {

    take:
    ch_raw_fastq_gz // channel: [ [ meta.id, meta.run, meta.group ], [ reads_1.fastq.gz, reads_2.fastq.gz ] ]

    main:
    ch_versions             = Channel.empty()
    ch_multiqc_files        = Channel.empty()
    ch_workdirs_to_clean    = Channel.empty()

    //
    // MODULE: Run FastQC on raw reads
    //
    FASTQC_RAW (
        ch_raw_fastq_gz
    )
    ch_multiqc_files    = ch_multiqc_files.mix(FASTQC_RAW.out.zip.collect{it[1]})
    ch_versions         = ch_versions.mix(FASTQC_RAW.out.versions.first())

    /*----------------------------------------------------------------------------
        Read Preprocessing
    ------------------------------------------------------------------------------*/
    if (params.run_fastp) {
        //
        // MODULE: Run fastp on raw reads
        //
        FASTP(
            ch_raw_fastq_gz,
            [],
            false,
            false,
            false
        )
        ch_fastp_prefilt_fastq_gz   = FASTP.out.reads
        ch_versions                 = ch_versions.mix(FASTP.out.versions)
        ch_multiqc_files            = ch_multiqc_files.mix(FASTP.out.json.collect{ it[1] })

        // REMOVE EMPTY FASTQ FILES FROM CHANNEL
        ch_fastp_fastq_gz = rmEmptyFastQs(ch_fastp_prefilt_fastq_gz)

        // IDENTIFY WORKDIRS TO CLEAN
        ch_pre_fastp_workdirs = getWorkDirs(
            ch_raw_fastq_gz,
            ch_fastp_fastq_gz
        )
        ch_workdirs_to_clean = ch_workdirs_to_clean.mix(ch_pre_fastp_workdirs)
    } else {
        ch_fastp_fastq_gz = ch_raw_fastq_gz
    }

    /*----------------------------------------------------------------------------
        Read merging
    ------------------------------------------------------------------------------*/
    if (params.perform_run_merging) {
        // prepare reads for concatenating within runs
        ch_reads_forcat = ch_fastp_fastq_gz
            .map { meta, reads -> [ meta - meta.subMap('run'), reads ] }
            .groupTuple()
            .branch {
                meta, reads ->
                    cat:      reads.size() >= 2 // SE: [ [ meta ], [ S1_R1, S2_R1 ] ]; PE: [ [ meta ], [ [ S1_R1, S1_R2 ], [ S2_R1, S2_R2 ] ] ]
                    skip_cat: true              // Can skip merging if only single lanes
            }

        //
        // MODULE: Concatenate reads across runs, within samples
        //
        CAT_RUNMERGE(
            ch_reads_forcat.cat.map { meta, reads -> [ meta, reads.flatten() ] }
        )
        ch_cat_reads_fastq_gz   = CAT_RUNMERGE.out.reads
        ch_versions             = ch_versions.mix(CAT_RUNMERGE.out.versions)

        // Ensure we don't have nests of nests so that structure is in form expected for assembly
        ch_reads_forcat_skipped = ch_reads_forcat.skip_cat
            .map { meta, reads ->
                def new_reads = meta.single_end ? reads[0] : reads.flatten()
                [ meta, new_reads ]
            }

        // Combine single run and multi-run-merged data
        ch_runmerged_fastq_gz = ch_cat_reads_fastq_gz.mix(ch_reads_forcat_skipped)

        // IDENTIFY WORKDIRS TO CLEAN
        ch_pre_merge_workdirs = getWorkDirs(
            ch_fastp_fastq_gz,
            ch_runmerged_fastq_gz
        )
        ch_workdirs_to_clean = ch_workdirs_to_clean.mix(ch_pre_merge_workdirs)
    } else {
        ch_runmerged_fastq_gz = ch_fastp_fastq_gz
    }

    /*----------------------------------------------------------------------------
        Host read removal
    ------------------------------------------------------------------------------*/
    if (params.run_bowtie2_host_removal) {
        // prepare host fasta and bowtie2 index
        if (params.bowtie2_igenomes_host) {
            ch_bowtie2_host_fasta   = Channel.value(
                [ [ id:'bowtie2_fasta' ], file(params.genomes[params.bowtie2_igenomes_host].fasta, checkIfExists: true) ]
            )
            ch_bowtie2_host_index   = Channel.value(
                [ [ id:'bowtie2_index' ], file(params.genomes[params.bowtie2_igenomes_host].bowtie2, checkIfExists: true) ]
            )
        } else {
            ch_bowtie2_host_fasta   = Channel.value(
                [ [ id:'bowtie2_fasta' ], file(params.bowtie2_custom_host_fasta, checkIfExists: true) ]
            )
            ch_bowtie2_host_index   = null
        }

        //
        // SUBWORKFLOW: Remove host reads using Bowtie2
        //
        FASTQ_HOSTREMOVAL_BOWTIE2(
            ch_fastp_fastq_gz,
            ch_bowtie2_host_fasta,
            ch_bowtie2_host_index
        )
        ch_bt2_prefilt_fastq_gz = FASTQ_HOSTREMOVAL_BOWTIE2.out.fastq_gz
        ch_versions             = ch_versions.mix(FASTQ_HOSTREMOVAL_BOWTIE2.out.versions)
        ch_multiqc_files        = ch_multiqc_files.mix(FASTQ_HOSTREMOVAL_BOWTIE2.out.mqc.collect{ it[1] })

        // REMOVE EMPTY FASTQ FILES FROM CHANNEL
        ch_bt2_fastq_gz = rmEmptyFastQs(ch_bt2_prefilt_fastq_gz)

        // IDENTIFY WORKDIRS TO CLEAN
        ch_pre_bt2_workdirs = getWorkDirs(
            ch_fastp_fastq_gz,
            ch_bt2_fastq_gz
        )
        ch_workdirs_to_clean = ch_workdirs_to_clean.mix(ch_pre_bt2_workdirs)
    } else {
        ch_bt2_fastq_gz  = ch_fastp_fastq_gz
    }

    /*----------------------------------------------------------------------------
        Preprocessing analysis
    ------------------------------------------------------------------------------*/
    if (params.perform_run_merging || params.run_fastp || params.run_bowtie2_host_removal) {
        //
        // MODULE: Run FastQC on preprocessed reads
        //
        FASTQC_PREPROCESSED (
            ch_bt2_fastq_gz
        )
        ch_multiqc_files    = ch_multiqc_files.mix(FASTQC_PREPROCESSED.out.zip.collect{it[1]})
        ch_versions         = ch_versions.mix(FASTQC_PREPROCESSED.out.versions.first())
    }

    /*----------------------------------------------------------------------------
        Estimate virus enrichment
    ------------------------------------------------------------------------------*/
    if (params.run_viromeqc) {
        // create channel from params.viromeqc_db
        if (!params.viromeqc_db) {
            ch_viromeqc_db  = null
        } else {
            ch_viromeqc_db  = Channel.value(
                file(params.viromeqc_db, checkIfExists:true)
            )
        }

        //
        // SUBWORKFLOW: Estimate virus enrichment with ViromeQC
        //
        FASTQ_VIRUSENRICHMENT_VIROMEQC(
            ch_bt2_fastq_gz,
            ch_viromeqc_db
        )
        ch_vqc_enrich_tsv   = FASTQ_VIRUSENRICHMENT_VIROMEQC.out.enrichment_tsv
        ch_versions         = ch_versions.mix(FASTQ_VIRUSENRICHMENT_VIROMEQC.out.versions)
    } else {
        ch_vqc_enrich_tsv   = []
    }

    emit:
    preprocessed_fastq_gz   = ch_bt2_fastq_gz   // channel: [ [ meta.id, meta.group ], [ reads_1.fastq.gz, reads_2.fastq.gz ] ]
    multiqc_files           = ch_multiqc_files  // channel: /path.to/multiqc_files
    versions                = ch_versions       // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
