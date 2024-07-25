#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HoffLab/phageannotator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/HoffLab/phageannotator
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// PLUGINS: Installed from nextflow plugins
//
include { paramsSummaryMap  } from 'plugin/nf-validation'

//
// MODULES
//
include { LOGAN_DOWNLOAD    } from '../modules/local/logan/download/main'
include { MULTIQC           } from './modules/nf-core/multiqc/main'

//
// SUBWORKFLOWS
//
include { CSV_SRADOWNLOAD_FETCHNGS  } from '../subworkflows/local/csv_sradownload_fetchngs/main'
include { PIPELINE_COMPLETION       } from './subworkflows/local/utils_vtdb_createdb_pipeline'
include { PIPELINE_INITIALISATION   } from './subworkflows/local/utils_vtdb_createdb_pipeline'
include { methodsDescriptionText    } from './subworkflows/local/utils_vtdb_createdb_pipeline'
include { paramsSummaryMultiqc      } from './subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML    } from './subworkflows/nf-core/utils_nfcore_pipeline'

//
// WORKFLOWS
//
include { FASTQ_READPREPROCESSING_FASTQ } from './workflows/fastq_readpreprocessing_fastq/main'
include { FASTQ_READASSEMBLY_FASTA      } from './workflows/fastq_readassembly_fasta/main'
include { FASTA_VIRUSMINE_FASTA         } from './workflows/fasta_virusmine_fasta/main'
include { FASTA_VIRUSANNOTATE_TSV       } from './workflows/fasta_virusannotate_tsv/main'
include { FASTQFASTA_VIRUSANALYZE_TSV   } from './workflows/fastqfasta_virusanalyze_tsv/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow HOFFLAB_PHAGEANNOTATOR {

    take:
    input_fastq_gz  // channel: [ [ meta.id, meta.run, meta.group ], [ reads_1.fastq.gz, reads_2.fastq.gz ] ]
    input_fasta_gz  // channel: [ [ meta.id, meta.run, meta.group ], fasta.gz ]

    main:
    ch_versions             = Channel.empty()
    ch_multiqc_files        = Channel.empty()

    /*----------------------------------------------------------------------------
        Load and parse accessions file
    ------------------------------------------------------------------------------*/
    //
    // Create channels from input file provided through params.input and params.assembly_input
    //
    ch_sra_accessions = Channel
        .fromSamplesheet("sra_accessions")

    //
    // SUBWORKFLOW: Download SRA reads
    //
    if (params.sra_accessions) {
        ch_sra_fastq_gz     = CSV_SRADOWNLOAD_FETCHNGS(ch_sra_accessions, params.sra_download_method).fastq
        ch_versions         = ch_versions.mix(CSV_SRADOWNLOAD_FETCHNGS.out.versions)
        ch_input_fastq_gz   = input_fastq_gz.mix(ch_sra_fastq_gz)
    } else {
        ch_input_fastq_gz   = input_fastq_gz
    }


    //
    // WORKFLOW: Run read preprocessing
    //
    FASTQ_READPREPROCESSING_FASTQ(
        ch_input_fastq_gz
    )
    ch_preprocessed_fastq_gz = FASTQ_READPREPROCESSING_FASTQ.out.preprocessed_fastq_gz
    ch_versions              = ch_versions.mix(FASTQ_READPREPROCESSING_FASTQ.out.versions)
    ch_multiqc_files         = ch_multiqc_files.mix(FASTQ_READPREPROCESSING_FASTQ.out.multiqc_files)


    //
    // MODULE: Download Logan assemblies
    //
    if (params.logan_accessions) {
        // load logan accessions from input samplesheet
        ch_logan_accession_files = Channel
            .fromFilePairs(params.logan_accessions, size: 1)
            .map { meta, accession_file ->
                [ [ id:"logan_acc_" + meta, assembler:"logan" ], accession_file ]
            }

        //
        // MODULE: Download Logan assemblies
        //
        LOGAN_DOWNLOAD(
            ch_logan_accession_files
        )
        ch_logan_fasta_zst  = LOGAN_DOWNLOAD.out.fasta
        ch_input_fasta_gz   = input_fasta_gz.mix(ch_logan_fasta_zst)
        ch_versions         = ch_versions.mix(LOGAN_DOWNLOAD.out.versions)
    } else {
        ch_input_fasta_gz   = input_fasta_gz
    }


    //
    // WORKFLOW: Run read assembly
    //
    FASTQ_READASSEMBLY_FASTA(
        ch_preprocessed_fastq_gz,
        ch_input_fasta_gz
    )
    ch_assemblies_fasta_gz  = FASTQ_READASSEMBLY_FASTA.out.assemblies_fasta_gz
    ch_versions             = ch_versions.mix(FASTQ_READASSEMBLY_FASTA.out.versions)
    ch_multiqc_files        = ch_multiqc_files.mix(FASTQ_READASSEMBLY_FASTA.out.multiqc_files)


    //
    // WORKFLOW: Mine viral sequences
    //
    FASTA_VIRUSMINE_FASTA(
        ch_preprocessed_fastq_gz,
        ch_assemblies_fasta_gz
    )
    ch_viruses_fasta_gz = FASTA_VIRUSMINE_FASTA.out.viruses_fasta_gz
    ch_versions         = ch_versions.mix(FASTA_VIRUSMINE_FASTA.out.versions)
    ch_multiqc_files    = ch_multiqc_files.mix(FASTA_VIRUSMINE_FASTA.out.multiqc_files)



    //
    // WORKFLOW: Cluster viral sequences
    //
    VIRUS_CLUSTER(
        ch_filt_viruses_fasta_gz
    )
    ch_cluster_repts_fasta_gz   = VIRUS_CLUSTER.out.filt_viruses_fasta_gz
    ch_versions                 = ch_versions.mix(VIRUS_CLUSTER.out.versions)
    ch_multiqc_files            = ch_multiqc_files.mix(VIRUS_CLUSTER.out.multiqc_files)


    //
    // WORKFLOW: Annotate viral sequences
    //
    VIRUS_ANNOTATE(
        ch_filt_viruses_fasta_gz
    )
    ch_cluster_repts_fasta_gz   = VIRUS_ANNOTATE.out.filt_viruses_fasta_gz
    ch_versions                 = ch_versions.mix(VIRUS_ANNOTATE.out.versions)
    ch_multiqc_files            = ch_multiqc_files.mix(VIRUS_ANNOTATE.out.multiqc_files)


    //
    // WORKFLOW: Analyze viral sequences
    //
    READVIRUS_ANALYZE(
        ch_filt_viruses_fasta_gz
    )
    ch_cluster_repts_fasta_gz   = READVIRUS_ANALYZE.out.filt_viruses_fasta_gz
    ch_versions                 = ch_versions.mix(READVIRUS_ANALYZE.out.versions)
    ch_multiqc_files            = ch_multiqc_files.mix(READVIRUS_ANALYZE.out.multiqc_files)


    /*----------------------------------------------------------------------------
        Report generation
    ------------------------------------------------------------------------------*/
    // Collate and save software versions
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    // Prepare MultiQC inputs
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

    //
    // MODULE: MultiQC
    //
    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    ch_multiqc_report  = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html

    emit:
    multiqc_report = ch_multiqc_report

}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    HOFFLAB_PHAGEANNOTATOR (
        PIPELINE_INITIALISATION.out.fastqs,
        PIPELINE_INITIALISATION.out.fastas
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        HOFFLAB_PHAGEANNOTATOR.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
