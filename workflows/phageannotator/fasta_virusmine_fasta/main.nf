//
// WORKFLOW: Identify and QC viral sequences in assemblies
//


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// FUNCTIONS: Local functions
//
include { getWorkDirs   } from '../../lib/NfCleanUp.groovy'

//
// MODULES: Local modules
//
include { TANTAN    } from '../../modules/local/tantan/main'

//
// SUBWORKFLOWS: Consisting of a mix of local and nf-core/modules
//
include { FASTA_VIRUSCLASSIFICATION_GENOMAD } from '../../subworkflows/local/fasta_virusclassification_genomad/main'
include { FASTQFASTA_VIRUSEXTENSION_COBRA   } from '../../subworkflows/local/fastqfasta_virusextension_cobra/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// PLUGINS
//

//
// MODULES: Installed directly from nf-core/modules
//
include { GENOMAD_ENDTOEND as GENOMAD_COBRA     } from '../../modules/nf-core/genomad/endtoend/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow FASTA_VIRUSMINE_FASTA {


    take:
    ch_preprocessed_fastq_gz    // channel: [ [ meta.id, meta.group ], [ reads_1.fastq.gz, reads_2.fastq.gz ] ]
    ch_assemblies_fasta_gz      // channel: [ [ meta.id, meta.group, meta.assembler ], fasta.gz ]


    main:
    ch_versions             = Channel.empty()
    ch_multiqc_files        = Channel.empty()
    ch_workdirs_to_clean    = Channel.empty()


    /*----------------------------------------------------------------------------
        Identify DTR/ITRs
    ------------------------------------------------------------------------------*/
    if (params.run_trfinder) {
        ch_trfinder_input = ch_assemblies_fasta_gz
            .map { meta, fasta ->
                meta.assembler = meta.assembler + "_trfinder"
                [ meta, fasta ]
            }
        //
        // MODULE: Find terminal repeats in fasta files
        //
        TRFINDER(
            ch_trfinder_input
        )
        ch_trfinder_fasta_gz    = TRFINDER.out.tr_fasta
        ch_trfinder_stats_tsv   = TRFINDER.out.tr_stats
        ch_versions             = ch_versions.mix(TRFINDER.out.versions)

        // delete pre trfinder fasta to save space
        ch_pre_trfinder_workdirs    = getWorkDirs(
            ch_assemblies_fasta_gz,
            ch_trfinder_fasta_gz,
            []
        )
        ch_workdirs_to_clean        = ch_workdirs_to_clean.mix(ch_pre_trfinder_workdirs)
    } else {
        ch_trfinder_fasta_gz    = ch_assemblies_fasta_gz
        ch_trfinder_stats_tsv   = ch_assemblies_fasta_gz
            .map { meta, fasta -> [ meta, [] ] }
    }


    /*----------------------------------------------------------------------------
        De novo virus classification
    ------------------------------------------------------------------------------*/
    if (params.run_genomad || params.run_cobra) {
        // create channel from params.genomad_db
        if (!params.genomad_db){
            ch_genomad_db   = null
        } else {
            ch_genomad_db   = Channel.value(
                file(params.genomad_db, checkIfExists:true)
            )
        }

        //
        // SUBWORKFLOW: Download and run geNomad
        //
        FASTA_VIRUSCLASSIFICATION_GENOMAD(
            ch_trfinder_fasta_gz,
            ch_genomad_db
        )
        ch_genomad_summary_tsv      = FASTA_VIRUSCLASSIFICATION_GENOMAD.out.virus_summary_tsv
        ch_genomad_scores_tsv       = FASTA_VIRUSCLASSIFICATION_GENOMAD.out.agg_class
        ch_genomad_taxonomy_tsv     = FASTA_VIRUSCLASSIFICATION_GENOMAD.out.taxonomy
        ch_genomad_genes_tsv        = FASTA_VIRUSCLASSIFICATION_GENOMAD.out.genes
        ch_genomad_prefilt_fasta_gz = FASTA_VIRUSCLASSIFICATION_GENOMAD.out.virus_fasta_gz
        ch_versions                 = ch_versions.mix(FASTA_VIRUSCLASSIFICATION_GENOMAD.out.versions)

        // remove empty fastA files from channel
        if (params.use_genomad_fasta) {
            ch_genomad_fasta_gz = ch_genomad_prefilt_fasta_gz
                .filter { meta, fasta ->
                    try {
                        fasta.countFasta(limit: 5) > 1
                    } catch (EOFException) {
                        log.warn "[HoffLab/phageannotator]: ${fasta} has an EOFException, this is likely an empty gzipped file."
                    }
                }
            if (!params.run_cobra) {
                // identify intermediate workDirs to clean
                ch_pre_genomad_workdirs = getWorkDirs(
                    ch_seqkit_seq_fasta_gz,
                    ch_genomad_fasta_gz,
                    []
                )
                ch_workdirs_to_clean    = ch_workdirs_to_clean.mix(ch_pre_genomad_workdirs)
            }
        } else {
            ch_genomad_fasta_gz = ch_trfinder_fasta_gz
        }
    } else {
        ch_genomad_fasta_gz     = ch_trfinder_fasta_gz
        ch_genomad_scores_tsv   = ch_trfinder_fasta_gz.map { meta, fasta -> [ meta, [] ] }
        ch_genomad_taxonomy_tsv = ch_genomad_scores_tsv
        ch_genomad_genes_tsv    = ch_genomad_scores_tsv
    }


    /*----------------------------------------------------------------------------
        Extend viral contigs
    ------------------------------------------------------------------------------*/
    if (params.run_cobra) {
        // filter to only single assemblies with reads available
        ch_cobra_input = ch_genomad_fasta_gz
            .join(ch_assembly_logs)
            .map { meta, fasta -> [ meta.id, meta, fasta ] }
            .combine( ch_bt2_fastq_gz.map { meta, fastq -> [ meta.id, meta, fastq ] }, by:0 )
            .map {
                id, meta_fasta, fasta, meta_fastq, fastq ->
                    meta_fasta.assembler =  meta_fasta.assember + '_cobra'
                [ meta_fasta, fastq, fasta ]
            }
            .filter { meta, fastq, fasta ->
                meta.assembler.contains('metaspades_single') || meta.assembler.contains('megahit_single')
            }
            .multiMap { meta, fastq, fasta ->
                fastq: [ meta, fastq ]
                fasta: [ meta, fasta ]
            }

        //
        // SUBWORKFLOW: Extend assembled contigs
        //
        FASTQFASTA_VIRUSEXTENSION_COBRA(
            ch_cobra_input.fastq,
            ch_cobra_input.fasta,
            FASTA_VIRUSCLASSIFICATION_GENOMAD.out.virus_summary_tsv
        )
        ch_cobra_prefilt_fasta_gz   = FASTQFASTA_VIRUSEXTENSION_COBRA.out.extended_fasta
        ch_cobra_summary_tsv        = FASTQFASTA_VIRUSEXTENSION_COBRA.out.cobra_summary_tsv
        ch_versions                 = ch_versions.mix(FASTQFASTA_VIRUSEXTENSION_COBRA.out.versions)

        // remove empty fastA files from channel
        ch_cobra_filt_fasta_gz = ch_cobra_prefilt_fasta_gz
            .filter { meta, fasta ->
                try {
                    fasta.countFasta(limit: 5) > 1
                } catch (EOFException) {
                    log.warn "[HoffLab/phageannotator]: ${fasta} has an EOFException, this is likely an empty gzipped file."
                }
            }

        //
        // MODULE: Run geNomad on extended viral contigs
        //
        GENOMAD_COBRA(
            ch_cobra_filt_fasta_gz,
            FASTA_VIRUSCLASSIFICATION_GENOMAD.out.genomad_db.first()
        )
        ch_cobra_genomad_fasta_gz = GENOMAD_COBRA.out.virus_fasta

        ch_genomad_summary_tsv  = ch_genomad_summary_tsv.mix(GENOMAD_COBRA.out.virus_summary)
        ch_genomad_scores_tsv   = GENOMAD_COBRA.out.agg_class
        ch_versions             = ch_versions.mix(GENOMAD_COBRA.out.versions)

        // combine cobra extended viruses with unextended viruses
        ch_cobra_fasta_gz = ch_genomad_fasta_gz.mix(ch_cobra_genomad_fasta_gz)
    } else {
        ch_cobra_fasta_gz       = ch_genomad_fasta_gz
        ch_cobra_summary_tsv    = Channel.empty()
    }


    /*----------------------------------------------------------------------------
        Assess virus quality and filter
    ------------------------------------------------------------------------------*/
    if (params.run_checkv || params.run_nucleotide_stats) {
        // create channel from params.checkv_db
        if (!params.checkv_db){
            ch_checkv_db    = null
        } else {
            ch_checkv_db    = Channel.value(
                file(params.checkv_db, checkIfExists:true)
            )
        }

        //
        // SUBWORKFLOW: Assess virus quality with CheckV
        //
        FASTA_VIRUSQUALITY_CHECKV(
            ch_viruses_fasta_gz,
            ch_checkv_db
        )
        ch_quality_summary_tsv  = FASTA_VIRUSQUALITY_CHECKV.out.quality_summary
        ch_checkv_faa_gz        = FASTA_VIRUSQUALITY_CHECKV.out.proteins_faa_gz
        ch_checkv_fna_gz        = FASTA_VIRUSQUALITY_CHECKV.out.viruses_fna_gz
        ch_versions             = ch_versions.mix(FASTA_VIRUSQUALITY_CHECKV.out.versions)

        if (params.use_checkv_fasta) {
            ch_checkv_fasta_gz = ch_checkv_fna_gz

            // identify intermediate workDirs to clean
            ch_pre_checkv_workdirs  = getWorkDirs(
                ch_cobra_fasta_gz,
                ch_quality_summary_tsv,
                []
            )
            ch_workdirs_to_clean    = ch_workdirs_to_clean.mix(ch_pre_checkv_workdirs)
        } else {
            ch_checkv_fasta_gz = ch_viruses_fasta_gz
        }

    } else {
        ch_quality_summary_tsv  = []
        ch_checkv_fna_gz        = ch_viruses_fasta_gz
    }

    if (params.run_tantan) {
        //
        // MODULE: Identify low-complexity regions with tantan
        //
        TANTAN(
            ch_checkv_fna_gz
        )
        ch_tantan_tsv   = TANTAN.out.tantan
        ch_versions     = ch_versions.mix(TANTAN.out.versions)

        // identify intermediate workDirs to clean
        ch_pre_tantan_workdirs  = getWorkDirs(
            ch_checkv_fna_gz,
            ch_tantan_tsv,
            []
        )
        ch_workdirs_to_clean    = ch_workdirs_to_clean.mix(ch_pre_tantan_workdirs)
    } else {
        ch_tantan_tsv   = []
    }

    if (params.run_nucleotide_stats) {
        // join fasta and proteins for nucleotide stats input
        ch_nuc_stats_input  = ch_checkv_fna_gz
            .join(ch_checkv_faa_gz)
            .multiMap { it ->
                fasta: [ it[0], it[1] ]
                proteins: [ it[0], it[2] ]
            }

        //
        // MODULE: Calculate nucleotide stats
        //
        NUCLEOTIDESTATS(
            ch_nuc_stats_input.fasta,
            ch_nuc_stats_input.proteins
        )
        ch_nuc_stats_tsv    = NUCLEOTIDESTATS.out.nuc_stats
        ch_versions         = ch_versions.mix(NUCLEOTIDESTATS.out.versions)

        // identify intermediate workDirs to clean
        ch_pre_nuc_stats_workdirs   = getWorkDirs(
            ch_cobra_fasta_gz,
            ch_nuc_stats_tsv,
            []
        )
        ch_workdirs_to_clean        = ch_workdirs_to_clean.mix (ch_pre_nuc_stats_workdirs)
    } else {
        ch_nuc_stats_tsv    = []
    }


    /*----------------------------------------------------------------------------
        Filter sequences based on composition, completeness, and quality
    ------------------------------------------------------------------------------*/
    if (params.run_sequence_filtering) {
        // join all virus classification files for input
        ch_filt_input = ch_checkv_fna_gz
            .join(ch_genomad_summary_tsv)
            .join(ch_quality_summary_tsv)
            .join(ch_tantan_tsv)
            .join(ch_nuc_stats_tsv)
            .multiMap { it ->
                fasta:              [ it[0], it[1] ]
                genomad_summary:    [ it[0], it[2] ]
                quality_summary:    [ it[0], it[3] ]
                tantan:             [ it[0], it[4] ]
                nuc_stats:          [ it[0], it[5] ]
            }

        if (params.contigs_to_keep) {
            ch_contigs_to_keep_tsv  = Channel.fromPath(
                file(params.contigs_to_keep, checkIfExists: true)
                )
        } else {
            ch_contigs_to_keep_tsv  = []
        }

        //
        // MODULE: Filter virus sequences based on classification metrics
        //
        FILTERSEQUENCES(
            ch_filt_input.fasta,
            ch_filt_input.genomad_summary,
            ch_filt_input.quality_summary,
            ch_filt_input.tantan,
            ch_filt_input.nuc_stats,
            ch_contigs_to_keep_tsv
        )
        ch_filt_prefilt_fasta_gz    = FILTERSEQUENCES.out.fasta
        ch_filt_data_tsv            = FILTERSEQUENCES.out.tsv

        // remove empty fastA files from channel
        ch_filt_fasta_gz = ch_filt_prefilt_fasta_gz
        .filter { meta, fasta ->
            try {
                fasta.countFasta(limit: 5) > 1
            } catch (EOFException) {
                log.warn "[HoffLab/phageannotator]: ${fasta} has an EOFException, this is likely an empty gzipped file."
            }
        }

        // identify intermediate workDirs to clean
        ch_pre_filt_workdirs  = getWorkDirs(
            ch_checkv_fna_gz,
            ch_filt_fasta_gz,
            ".*_tantan.tsv",
        )
        ch_workdirs_to_clean    = ch_workdirs_to_clean.mix(ch_pre_filt_workdirs)
    } else {
        ch_filt_fasta_gz        = ch_checkv_fna_gz
    }


    emit:
    viruses_fasta_gz        = ch_cobra_fasta_gz         // channel: [ [ meta.id, meta.group, meta.assembler ], assembly.fasta.gz ]
    genomad_scores_tsv      = ch_genomad_scores_tsv     // channel: [ [ meta.id, meta.group ], [ scores.tsv ] ]
    genomad_taxonomy_tsv    = ch_genomad_taxonomy_tsv   // channel: [ [ meta.id, meta.group ], [ taxonomy.tsv ] ]
    genomad_genes_tsv       = ch_genomad_genes_tsv      // channel: [ [ meta.id, meta.group ], [ genes.tsv ] ]
    multiqc_files           = ch_multiqc_files          // channel: /path.to/multiqc_files
    versions                = ch_versions               // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
