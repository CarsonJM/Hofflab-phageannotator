/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// FUNCTIONS: Local functions
//
include { getWorkDirs                               } from '../../lib/NfCleanUp.groovy'

//
// MODULES: Local modules
//
// include { CLEANWORKDIRS                             } from '../../modules/local/cleanup/cleanworkdirs/main'
include { FILTERSEQUENCES                           } from '../../modules/local/filtersequences/main'
include { NUCLEOTIDESTATS                           } from '../../modules/local/nucleotidestats/main'
include { PHIST                                     } from '../../modules/local/phist/main'
include { SEQKIT_SEQ                                } from '../../modules/local/seqkit/seq/main'
// include { SKANI_TRIANGLE                            } from '../../modules/local/skani/triangle/main'
include { TANTAN                                    } from '../../modules/local/tantan/main'

//
// SUBWORKFLOWS: Consisting of a mix of local and nf-core/modules
//
// include { FASTA_PHAGEFUNCTION_PHAROKKA      } from '../../subworkflows/local/fasta_phagefunction_pharokka/main'
// include { FASTA_PHAGEHOST_IPHOP             } from '../../subworkflows/local/fasta_phagehost_iphop/main'
include { FASTA_VIRUSCLASSIFICATION_GENOMAD } from '../../subworkflows/local/fasta_virusclassification_genomad/main'
include { FASTA_VIRUSQUALITY_CHECKV         } from '../../subworkflows/local/fasta_virusquality_checkv/main'
include { FASTQ_HOSTREMOVAL_BOWTIE2         } from '../../subworkflows/local/fastq_hostremoval_bowtie2/main'
include { FASTQ_VIRUSENRICHMENT_VIROMEQC    } from '../../subworkflows/local/fastq_virusenrichment_viromeqc/main'
include { FASTQFASTA_VIRUSEXTENSION_COBRA   } from '../../subworkflows/local/fastqfasta_virusextension_cobra/main'
include { methodsDescriptionText            } from '../../subworkflows/local/utils_nfcore_phageannotator_pipeline'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// PLUGINS
//
include { paramsSummaryMap                      } from 'plugin/nf-validation'

//
// MODULES: Installed directly from nf-core/modules
//
// include { BACPHLIP                              } from '../../modules/nf-core/bacphlip/main'
include { CAT_FASTQ as CAT_RUNMERGE         } from '../../modules/nf-core/cat/fastq/main'
include { CAT_FASTQ as CAT_COASSEMBLY       } from '../../modules/nf-core/cat/fastq/main'
include { FASTP                             } from '../../modules/nf-core/fastp/main'
include { FASTQC as FASTQC_RAW              } from '../../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_PREPROCESSED     } from '../../modules/nf-core/fastqc/main'
include { GENOMAD_ENDTOEND as GENOMAD_COBRA } from '../../modules/nf-core/genomad/endtoend/main'
include { MULTIQC                           } from '../../modules/nf-core/multiqc/main'
include { SPADES as METASPADES_COASSEMBLY   } from '../../modules/nf-core/spades/main'
include { SPADES as METASPADES_SINGLE       } from '../../modules/nf-core/spades/main'

//
// SUBWORKFLOWS: Installed directly from nf-core/modules
//
include { paramsSummaryMultiqc                  } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                } from '../../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PHAGEANNOTATOR {

    take:
    input_reads_fastq_gz        // channel: [ [ meta.id, meta.group ], [ reads_1.fastq.gz, reads_2.fastq.gz ] ]
    input_assemblies_fasta_gz   // channel: [ [ meta.id, meta.group ], assembly.fasta ]

    main:
    ch_versions             = Channel.empty()
    ch_multiqc_files        = Channel.empty()
    ch_workdirs_to_clean    = Channel.empty()

    //
    // MODULE: Run FastQC on raw reads
    //
    FASTQC_RAW (
        input_reads_fastq_gz
    )
    ch_multiqc_files    = ch_multiqc_files.mix(FASTQC_RAW.out.zip.collect{it[1]})
    ch_versions         = ch_versions.mix(FASTQC_RAW.out.versions.first())


    /*----------------------------------------------------------------------------
        Read merging
    ------------------------------------------------------------------------------*/
    if (params.perform_run_merging) {
        // prepare reads for concatenating within runs
        ch_reads_forcat             = input_reads_fastq_gz
            .map {
                meta, reads ->
                    def meta_new    = meta - meta.subMap('run')
                [ meta_new, reads ]
            }
            .groupTuple()
            .branch {
                meta, reads ->
                    cat:      reads.size() >= 2 // SE: [ [ meta ], [ S1_R1, S2_R1 ] ]; PE: [ [ meta ], [ [ S1_R1, S1_R2 ], [ S2_R1, S2_R2 ] ] ]
                    skip_cat: true              // Can skip merging if only single lanes
            }

        ch_cat_reads_fastq_gz       = CAT_RUNMERGE(
            ch_reads_forcat.cat.map { meta, reads -> [ meta, reads.flatten() ] }
        ).reads

        // Ensure we don't have nests of nests so that structure is in form expected for assembly
        ch_reads_forcat_skipped     = ch_reads_forcat.skip_cat
            .map { meta, reads ->
                def new_reads = meta.single_end ? reads[0] : reads.flatten()
                [ meta, new_reads ]
            }

        // Combine single run and multi-run-merged data
        ch_merged_reads_fastq_gz    = ch_cat_reads_fastq_gz.mix(ch_reads_forcat_skipped)
        ch_versions                 = ch_versions.mix(CAT_RUNMERGE.out.versions)

        // identify workDirs to clean
        ch_pre_merge_workdirs       = getWorkDirs(
            input_reads_fastq_gz,
            ch_merged_reads_fastq_gz,
            []
        )
        ch_workdirs_to_clean        = ch_workdirs_to_clean.mix(ch_pre_merge_workdirs)
    } else {
        ch_merged_reads_fastq_gz = input_reads_fastq_gz
    }


    /*----------------------------------------------------------------------------
        Read Preprocessing
    ------------------------------------------------------------------------------*/
    if (params.run_fastp) {
        //
        // MODULE: Run fastp on raw reads
        //
        ch_fastp_fastq_gz = FASTP(
            ch_merged_reads_fastq_gz,
            [],
            false,
            false
        ).reads
        ch_versions     = ch_versions.mix(FASTP.out.versions)

        // identify workDirs to clean
        ch_pre_fastp_workdirs   = getWorkDirs(
            ch_merged_reads_fastq_gz,
            ch_fastp_fastq_gz,
            []
        )
        ch_workdirs_to_clean    = ch_workdirs_to_clean.mix(ch_pre_fastp_workdirs)
    } else {
        ch_fastp_fastq_gz = ch_merged_reads_fastq_gz
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
        ch_bt2_fastq_gz = FASTQ_HOSTREMOVAL_BOWTIE2(
            ch_fastp_fastq_gz,
            ch_bowtie2_host_fasta,
            ch_bowtie2_host_index
        ).fastq_gz
        ch_versions     = ch_versions.mix(FASTQ_HOSTREMOVAL_BOWTIE2.out.versions)

        // identify workDirs to clean
        ch_pre_bt2_workdirs     = getWorkDirs(
            ch_fastp_fastq_gz,
            ch_bt2_fastq_gz,
            []
        )
        ch_workdirs_to_clean    = ch_workdirs_to_clean.mix(ch_pre_bt2_workdirs)
    } else {
        ch_bt2_fastq_gz  = ch_fastp_fastq_gz
    }


    /*----------------------------------------------------------------------------
        Estimate virus enrichment
    ------------------------------------------------------------------------------*/
    if (params.run_viromeqc) {
        // create channel from params.viromeqc_db
        if (!params.viromeqc_db){
            ch_viromeqc_db  = null
        } else {
            ch_viromeqc_db  = Channel.value(
                file(params.viromeqc_db, checkIfExists:true)
            )
        }

        //
        // SUBWORKFLOW: Estimate virus enrichment with ViromeQC
        //
        ch_vqc_enrich_tsv   = FASTQ_VIRUSENRICHMENT_VIROMEQC(ch_bt2_fastq_gz, ch_viromeqc_db).enrichment_tsv
        ch_versions         = ch_versions.mix(FASTQ_VIRUSENRICHMENT_VIROMEQC.out.versions)
    } else {
        ch_vqc_enrich_tsv   = []
    }


    /*----------------------------------------------------------------------------
        Assemble reads into contigs/scaffolds
    ------------------------------------------------------------------------------*/
    if (params.run_metaspades_single) {
        // prepare reads for metaspades input
        ch_metaspades_single_input  = ch_bt2_fastq_gz.map { meta, fastq ->
            [ meta, fastq, [], [] ]
        }

        //
        // MODULE: Assemble reads with metaSPAdes
        //
        METASPADES_SINGLE(
            ch_metaspades_single_input,
            [],
            []
        )
        if (params.metaspades_use_scaffolds) {
            ch_metaspades_single_fasta_gz   = METASPADES_SINGLE.out.scaffolds.map { meta, fasta ->
                [ meta + [ coassembly: false ], fasta ]
            }
        } else {
            ch_metaspades_single_fasta_gz   = METASPADES_SINGLE.out.contigs.map { meta, fasta ->
                [ meta + [ coassembly: false ], fasta ]
            }
        }
        ch_metaspades_fasta_gz              = Channel.empty()
        ch_metaspades_fasta_gz              = ch_metaspades_fasta_gz.mix(ch_metaspades_single_fasta_gz)
        ch_versions                         = ch_versions.mix(METASPADES_SINGLE.out.versions)

        // identify workDirs to clean
        ch_pre_spades_workdirs  = getWorkDirs(
            ch_bt2_fastq_gz,
            ch_metaspades_single_fasta_gz,
            []
        )
        ch_workdirs_to_clean    = ch_workdirs_to_clean.mix(ch_pre_spades_workdirs)
    } else {
        ch_metaspades_fasta_gz  = Channel.empty()
    }

    if (params.run_metaspades_coassembly) {
        // group and set group as new id
        ch_cat_coassembly_fastq_gz      = ch_bt2_fastq_gz
            .map { meta, reads -> [ meta.group, meta, reads ] }
            .groupTuple(by: 0, sort:'deep')
            .map { group, metas, reads ->
                def meta                = [:]
                meta.id                 = "group-$group"
                meta.group              = group
                meta.single_end         = params.single_end
                if ( params.single_end ) [ meta, reads.collect { it }, [] ]
                else [ meta, reads.flatten() ]
            }

        //
        // MODULE: Combine reads within groups for coassembly
        //
        ch_coassembly_fastq_gz  = CAT_COASSEMBLY(
            ch_cat_coassembly_fastq_gz
        ).reads
        ch_versions             = ch_versions.mix(CAT_COASSEMBLY.out.versions)

        // prepare reads for metaspades input
        ch_metaspades_coassembly_input  = ch_coassembly_fastq_gz.map { meta, fastq ->
            [ meta, fastq, [], [] ]
        }

        //
        // MODULE: Assemble reads with metaSPAdes
        //
        METASPADES_COASSEMBLY(
            ch_metaspades_coassembly_input,
            [],
            []
        )
        if (params.metaspades_use_scaffolds) {
            ch_metaspades_co_fasta_gz   = METASPADES_COASSEMBLY.out.scaffolds.map { meta, fasta ->
                [ meta + [ coassembly: true ], fasta ]
            }
        } else {
            ch_metaspades_co_fasta_gz   = METASPADES_COASSEMBLY.out.contigs.map { meta, fasta ->
                [ meta + [ coassembly: true ], fasta ]
            }
        }
        ch_metaspades_fasta_gz          = ch_metaspades_fasta_gz.mix(ch_metaspades_co_fasta_gz)
        ch_versions                     = ch_versions.mix(METASPADES_COASSEMBLY.out.versions)

        // identify workDirs to clean
        ch_pre_spades_workdirs  = getWorkDirs(
            ch_bt2_fastq_gz,
            ch_metaspades_fasta_gz,
            []
        )
        ch_workdirs_to_clean    = ch_workdirs_to_clean.mix(ch_pre_spades_workdirs)
    }


    /*----------------------------------------------------------------------------
        Run quick quality analysis of assemblies
    ------------------------------------------------------------------------------*/
    if (params.run_seqkit_stats) {
        //
        // MODULE: Filter assemblies by length
        //
        ch_seqkit_stats_tsv = SEQKIT_STATS(ch_metaspades_fasta_gz).fastx
        ch_versions         = ch_versions.mix(SEQKIT_SEQ.out.versions)
    } else {
        ch_seqkit_stats_tsv = []
    }


    /*----------------------------------------------------------------------------
        Remove low-length assemblies
    ------------------------------------------------------------------------------*/
    if (params.run_seqkit_seq) {
        //
        // MODULE: Filter assemblies by length
        //
        ch_seqkit_seq_fasta_gz  = SEQKIT_SEQ(ch_metaspades_fasta_gz).fastx
        ch_versions             = ch_versions.mix(SEQKIT_SEQ.out.versions)

        if (!params.run_cobra) {
            // identify workDirs to clean
            ch_pre_seqkit_seq_workdirs  = getWorkDirs (
                ch_metaspades_fasta_gz,
                ch_seqkit_seq_fasta_gz,
                []
            )
            ch_workdirs_to_clean        = ch_workdirs_to_clean.mix(ch_pre_seqkit_seq_workdirs)
        }
    } else {
        ch_seqkit_seq_fasta_gz  = ch_metaspades_fasta_gz
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
        ch_genomad_summary_tsv      = FASTA_VIRUSCLASSIFICATION_GENOMAD(ch_seqkit_seq_fasta_gz, ch_genomad_db).virus_summary_tsv
        ch_genomad_fasta_gz         = FASTA_VIRUSCLASSIFICATION_GENOMAD.out.virus_fasta_gz
        ch_versions                 = ch_versions.mix(FASTA_VIRUSCLASSIFICATION_GENOMAD.out.versions)

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
        ch_genomad_fasta_gz     = ch_seqkit_seq_fasta_gz
        ch_genomad_summary_tsv  = []
    }


    /*----------------------------------------------------------------------------
        Extend viral contigs
    ------------------------------------------------------------------------------*/
    if (params.run_cobra) {
        // filter to only single assemblies with reads available
        ch_cobra_input = ch_bt2_fastq_gz
            .map { meta, fastq -> [ meta.id, meta, fastq ] }
            .join(ch_metaspades_fasta_gz.map { meta, fasta -> [ meta.id, meta, fasta ] } )
            .map {
                id, meta_fastq, fastq, meta_fasta, fasta ->
                [ meta_fasta, fastq, fasta ]
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
            FASTA_VIRUSCLASSIFICATION_GENOMAD.out.virus_summary_tsv,
            "metaspades",
            params.cobra_mink,
            params.cobra_maxk
        )
        ch_extended_fasta_gz    = FASTQFASTA_VIRUSEXTENSION_COBRA.out.extended_fasta
            .map { meta, fasta ->
                [ meta + [ coassembly: false, extended: true ], fasta ]
            }
        ch_cobra_summary_tsv    = FASTQFASTA_VIRUSEXTENSION_COBRA.out.cobra_summary_tsv
        ch_versions             = ch_versions.mix(FASTQFASTA_VIRUSEXTENSION_COBRA.out.versions)

        //
        // SUBWORKFLOW: Run geNomad on extended viral contigs
        //
        GENOMAD_COBRA(ch_extended_fasta_gz, FASTA_VIRUSCLASSIFICATION_GENOMAD.out.genomad_db)
        ch_genomad_summary_tsv  = ch_genomad_summary_tsv.map { meta, summary ->
                [ meta + [ extended: false ], summary ]
            }
            .mix(GENOMAD_COBRA.out.virus_summary)
        ch_versions             = ch_versions.mix(GENOMAD_COBRA.out.versions)

        // replace single assemblies with cobra extended assemblies
        ch_cobra_fasta_gz = ch_genomad_fasta_gz
            .filter { meta, fasta -> meta.coassembly == true }
            .map { meta, fasta -> [ meta + [ extended: false ], fasta ] }
            .mix(GENOMAD_COBRA.out.virus_fasta)
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
        ch_quality_summary_tsv  = FASTA_VIRUSQUALITY_CHECKV(ch_cobra_fasta_gz, ch_checkv_db).quality_summary_tsv
        ch_checkv_faa_gz        = FASTA_VIRUSQUALITY_CHECKV.out.proteins_faa_gz
        ch_checkv_fna_gz        = FASTA_VIRUSQUALITY_CHECKV.out.viruses_fna_gz
        ch_versions             = ch_versions.mix(FASTA_VIRUSQUALITY_CHECKV.out.versions)

        // identify intermediate workDirs to clean
        ch_pre_checkv_workdirs  = getWorkDirs(
            ch_cobra_fasta_gz,
            ch_quality_summary_tsv,
            []
        )
        ch_workdirs_to_clean    = ch_workdirs_to_clean.mix(ch_pre_checkv_workdirs)
    } else {
        ch_quality_summary_tsv  = []
        ch_checkv_fna_gz        = ch_cobra_fasta_gz
    }

    if (params.run_tantan) {
        //
        // MODULE: Identify low-complexity regions with tantan
        //
        ch_tantan_tsv   = TANTAN(ch_checkv_fna_gz).tantan
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
        ch_nuc_stats_tsv    = NUCLEOTIDESTATS(ch_nuc_stats_input.fasta, ch_nuc_stats_input.proteins).nuc_stats
        ch_versions         = ch_versions.mix(NUCLEOTIDESTATS.out.versions)

        // identify intermediate workDirs to clean
        ch_pre_nuc_stats_workdirs   = getWorkDirs(
            ch_cobra_fasta_gz,
            ch_nuc_stats_tsv,
            []
        )
        ch_workdirs_to_clean        = ch_workdirs_to_clean.mix ( ch_pre_nuc_stats_workdirs )
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
        ch_filt_fasta_gz    = FILTERSEQUENCES.out.fasta
        ch_filt_data_tsv    = FILTERSEQUENCES.out.tsv

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


    /*----------------------------------------------------------------------------
        Predict phage hosts
    ------------------------------------------------------------------------------*/
    if (params.run_iphop){
        // create channel from params.checkv_db
        if (!params.iphop_db){
            ch_iphop_db = null
        } else {
            ch_iphop_db = file(params.iphop_db, checkIfExists:true)
        }

        //
        // SUBWORKFLOW: Download database and predict phage hosts
        //
        ch_iphop_predictions_tsv    = FASTA_PHAGE_HOST_IPHOP(ch_filt_fasta_gz, ch_iphop_db).host_predictions_tsv
        ch_versions                 = ch_versions.mix(FASTA_PHAGE_HOST_IPHOP.out.versions)
    } else {
        ch_iphop_predictions_tsv    = Channel.empty()
    }

    if (params.run_phist){
        // Split phages into individual fasta files for phist
        ch_phist_fastas = PHIST_SPLIT(ch_trfinder_fasta_gz).dir
        ch_versions     = ch_versions.mix(PHIST_SPLIT.out.versions)

        //
        // MODULE: Run PHIST to predict phage host using shared kmers
        //
        ch_phist_predictions_tsv    = PHIST(ch_phist_fastas, params.phist_bacteria_db_dir).predictions
        ch_phist_comm_kmers_tsv     = PHIST.out.common_kmers
        ch_versions                 = ch_versions.mix(PHIST.out.versions )
    } else {
        ch_iphop_predictions_tsv    = Channel.empty()
    }

    if (params.run_crispr_blast){
        // TODO: Add crispr blast option
    } else {
        ch_crispr_blast_tsv = Channel.empty()
    }


    /*----------------------------------------------------------------------------
        Phage functional annotation
    ------------------------------------------------------------------------------*/
    if (params.run_pharokka) {
        // create channel from params.pharokka_db
        if (!params.pharokka_db){
            ch_pharokka_db = null
        } else {
            ch_pharokka_db = Channel.value(
                file( params.pharokka_db, checkIfExists:true )
                )
        }

        //
        // SUBWORKFLOW: Functionally annotate phage sequences
        //
        ch_pharokka_gbk_gz      = FASTA_PHAGE_FUNCTION_PHAROKKA(ch_anicluster_reps_fasta, ch_pharokka_db).pharokka_gbk_gz
        ch_pharokka_output_tsv  = FASTA_PHAGE_FUNCTION_PHAROKKA.out.pharokka_final_output_tsv
        ch_versions             = ch_versions.mix(FASTA_PHAGE_FUNCTION_PHAROKKA.out.versions)
    } else {
        ch_pharokka_gbk_gz      = Channel.empty()
        ch_pharokka_output_tsv  = Channel.empty()
    }


    /*----------------------------------------------------------------------------
        Predict virus lifestyle
    ------------------------------------------------------------------------------*/
    if (params.run_bacphlip) {
        //
        // MODULE: Predict phage lifestyle with BACPHLIP
        //
        ch_bacphlip_lifestyle_tsv   = BACPHLIP(ch_anicluster_reps_fasta).bacphlip_results
        ch_versions                 = ch_versions.mix(BACPHLIP.out.versions)
    } else {
        ch_bacphlip_lifestyle_tsv   = Channel.empty()
    }

    // TODO: Add integration status check
    // TODO: Add PHROG integrase identification


    // /*----------------------------------------------------------------------------
    //     Predict if proviruses are active
    // ------------------------------------------------------------------------------*/
    // // TODO: Update python code with Adam's recommendations
    // if ( params.run_propagate ){
    //     ch_provirus_activity_tsv    = FASTQ_FASTA_PROVIRUS_ACTIVITY_PROPAGATE (
    //         ch_fastq_gz.fastq_included,
    //         fasta_gz,
    //         ch_virus_summaries_tsv,
    //         ch_quality_summary_tsv,
    //         ch_clusters_tsv,
    //         params.propagate_min_ani,
    //         params.propagate_min_qcov,
    //         params.propagate_min_tcov
    //         ).propagate_results_tsv
    //     ch_versions                 = ch_versions.mix ( FASTQ_FASTA_PROVIRUS_ACTIVITY_PROPAGATE.out.versions )
    // } else {
    //     ch_provirus_activity_tsv    = Channel.empty()
    // }


    // /*----------------------------------------------------------------------------
    //     Clean up intermediate files
    // ------------------------------------------------------------------------------*/
    // if (params.remove_intermediate_files) {
    //     //
    //     // MODULE: Clean up intermediate working directories
    //     //
    //     ch_workdirs_to_clean_unique = ch_workdirs_to_clean.unique()
    //     CLEANWORKDIRS (ch_workdirs_to_clean_unique)
    // }



    /*----------------------------------------------------------------------------
        Report generation
    ------------------------------------------------------------------------------*/
    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
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

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/