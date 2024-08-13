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
include { FASTG2GFA                             } from '../../modules/local/fastg2gfa/main'
include { PLASS_PENGUIN as PENGUIN_SINGLE       } from '../../modules/local/plass/penguin/main'
include { PLASS_PENGUIN as PENGUIN_COASSEMBLY   } from '../../modules/local/plass/penguin/main'
include { SEQKIT_SEQ                            } from '../../modules/local/seqkit/seq/main'
include { SEQKIT_SPLIT                          } from '../../modules/local/seqkit/split/main'
include { SEQKIT_STATS                          } from '../../modules/local/seqkit/stats/main'
include { TANTAN                                } from '../../modules/local/tantan/main'

//
// SUBWORKFLOWS: Consisting of a mix of local and nf-core/modules
//
include { FASTQGFA_VIRUSEXTENSION_PHABLES   } from '../../subworkflows/local/fastqgfa_virusextension_phables/main'
include { FASTQFASTA_VIRUSEXTENSION_COBRA   } from '../../subworkflows/local/fastqfasta_virusextension_cobra/main'

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
include { CAT_FASTQ as CAT_COASSEMBLY           } from '../../modules/nf-core/cat/fastq/main'
include { MEGAHIT as MEGAHIT_SINGLE             } from '../../modules/nf-core/megahit/main'
include { MEGAHIT as MEGAHIT_COASSEMBLY         } from '../../modules/nf-core/megahit/main'
include { SPADES as METASPADES_COASSEMBLY       } from '../../modules/nf-core/spades/main'
include { SPADES as METASPADES_SINGLE           } from '../../modules/nf-core/spades/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PHAGEANNOTATOR {

    take:
    ch_preprocessed_fastq_gz    // channel: [ [ meta.id, meta.group ], [ reads_1.fastq.gz, reads_2.fastq.gz ] ]
    ch_raw_fasta_gz             // channel: [ [ meta.id, meta.group ], fasta.gz ]

    main:
    ch_versions             = Channel.empty()
    ch_multiqc_files        = Channel.empty()
    ch_workdirs_to_clean    = Channel.empty()


    /*----------------------------------------------------------------------------
        Assemble reads into contigs/scaffolds
    ------------------------------------------------------------------------------*/
    if (
        params.run_metaspades_single ||
        params.run_megahit_single ||
        params.run_penguin_single ||
        params.run_metaspades_coassembly ||
        params.run_megahit_coassembly ||
        params.run_penguin_coassembly
    ) {
        ch_assemblies_prefilt_fasta_gz  = Channel.empty()
        ch_assembly_graph_gz            = Channel.empty()
        ch_assembly_log                 = Channel.empty()

        if (params.run_metaspades_single) {
            // prepare reads for metaspades input
            ch_metaspades_single_input = ch_preprocessed_fastq_gz
                .map { meta, fastq ->
                    [ meta + [ assembler: 'metaspades_single' ], fastq, [], [] ]
                }
                .branch { meta, fastq, extra1, extra2 ->
                    single_end: meta.single_end
                    paired_end: true
                }

            //
            // MODULE: Assemble reads individually with metaSPAdes
            //
            METASPADES_SINGLE(
                ch_metaspades_single_input.paired_end,
                [],
                []
            )
            if (params.metaspades_use_scaffolds) {
                ch_metaspades_single_fasta_gz   = METASPADES_SINGLE.out.scaffolds
            } else {
                ch_metaspades_single_fasta_gz   = METASPADES_SINGLE.out.contigs
            }
            ch_assemblies_prefilt_fasta_gz      = ch_assemblies_prefilt_fasta_gz.mix(ch_metaspades_single_fasta_gz)
            ch_assembly_graph_prefilt_gz        = ch_assembly_graph_gz.mix(METASPADES_SINGLE.out.gfa)
            ch_versions                         = ch_versions.mix(METASPADES_SINGLE.out.versions)
        }

        if (params.run_megahit_single) {
            //
            // MODULE: Assemble reads individually with MEGAHIT
            //
            MEGAHIT_SINGLE(
                ch_preprocessed_fastq_gz.map { meta, fasta -> [ meta + [ assembler: 'megahit_single' ], fasta ] }
            ).contigs
            ch_megahit_single_fasta_gz      = MEGAHIT_SINGLE.out.contigs
            ch_assemblies_prefilt_fasta_gz  = ch_assemblies_prefilt_fasta_gz.mix(ch_megahit_single_fasta_gz)
            ch_assembly_graph_prefilt_gz    = ch_assembly_graph_gz.mix(MEGAHIT_SINGLE.out.graph)
            ch_versions                     = ch_versions.mix(MEGAHIT_SINGLE.out.versions)
        }

        if (params.run_penguin_single) {
            //
            // MODULE: Assemble reads individually with PenguiN
            //
            PENGUIN_SINGLE(
                ch_preprocessed_fastq_gz.map { meta, fasta -> [ meta + [ assembler: 'penguin_single' ], fasta ] }
            )
            ch_penguin_single_fasta_gz      = PENGUIN_SINGLE.out.contigs
            ch_assemblies_prefilt_fasta_gz  = ch_assemblies_prefilt_fasta_gz.mix(ch_penguin_single_fasta_gz)
            ch_versions                     = ch_versions.mix(PENGUIN_SINGLE.out.versions)
        }

        if ( params.run_metaspades_coassembly || params.run_megahit_coassembly || params.run_penguin_coassembly ) {
            // group and set group as new id
            ch_cat_coassembly_fastq_gz = ch_bt2_fastq_gz
                .map { meta, reads -> [ meta.group, meta, reads ] }
                .groupTuple( by: 0, sort:'deep' )
                .map { group, meta, reads ->
                    def meta_new                = [:]
                    meta_new.id                 = "group-$group"
                    meta_new.group              = group
                    meta_new.single_end         = meta.single_end[0]
                    if ( meta_new.single_end ) {
                        return [ meta_new, reads.collect { it } ]
                    } else {
                        return [ meta_new, reads.flatten() ]
                    }
                }
                .branch {
                    meta, reads ->
                        coassembly: meta.single_end && reads.size() >= 2 || !meta.single_end && reads.size() >= 4
                        skip_coassembly: true   // Can skip coassembly if there is not multiple samples
                }
        }

        if ( params.run_metaspades_coassembly ) {
            //
            // MODULE: Combine reads within groups for coassembly
            //
            CAT_COASSEMBLY(
                ch_cat_coassembly_fastq_gz.coassembly
            )
            ch_coassembly_fastq_gz  = CAT_COASSEMBLY.out.reads
            ch_versions             = ch_versions.mix(CAT_COASSEMBLY.out.versions)

            // prepare reads for metaspades input
            ch_metaspades_coassembly_input = ch_coassembly_fastq_gz
                .map { meta, fastq ->
                    [ meta + [ assembler: 'metaspades_coassembly' ], fastq, [], [] ]
                }
                .branch { meta, fastq, extra1, extra2 ->
                    single_end: meta.single_end
                    paired_end: true
                }

            //
            // MODULE: Co-assemble reads with metaSPAdes
            //
            METASPADES_COASSEMBLY(
                ch_metaspades_coassembly_input.paired_end,
                [],
                []
            )
            if (params.metaspades_use_scaffolds) {
                ch_metaspades_co_fasta_gz   = METASPADES_COASSEMBLY.out.scaffolds
            } else {
                ch_metaspades_co_fasta_gz   = METASPADES_COASSEMBLY.out.contigs
            }
            ch_assemblies_prefilt_fasta_gz  = ch_assemblies_prefilt_fasta_gz.mix(ch_metaspades_co_fasta_gz)
            ch_assembly_graph_prefilt_gz    = ch_assembly_graph_gz.mix(METASPADES_COASSEMBLY.out.gfa)
            ch_versions                     = ch_versions.mix(METASPADES_COASSEMBLY.out.versions)
        }

        if (params.run_megahit_coassembly) {
            //
            // MODULE: Co-assemble reads with MEGAHIT
            //
            MEGAHIT_COASSEMBLY(
                ch_cat_coassembly_fastq_gz.coassembly.map{ meta, fasta -> [ meta + [ assembler: 'megahit_coassembly' ], fasta ] }
            )
            ch_megahit_co_fasta_gz          = MEGAHIT_COASSEMBLY.out.contigs
            ch_versions                     = ch_versions.mix(MEGAHIT_COASSEMBLY.out.versions)
            ch_assemblies_prefilt_fasta_gz  = ch_assemblies_prefilt_fasta_gz.mix(ch_megahit_co_fasta_gz)
            ch_assembly_graph_prefilt_gz    = ch_assembly_graph_gz.mix(MEGAHIT_COASSEMBLY.out.graph)
        }

        if (params.run_penguin_coassembly) {
            //
            // MODULE: Co-assemble reads with PenguiN
            //
            PENGUIN_COASSEMBLY(
                ch_cat_coassembly_fastq_gz.coassembly.map { meta, fasta -> [ meta + [ assembler: 'penguin_coassembly' ], fasta ] }
            )
            ch_penguin_co_fasta_gz          = PENGUIN_COASSEMBLY.out.contigs
            ch_versions                     = ch_versions.mix(PENGUIN_COASSEMBLY.out.versions)
            ch_assemblies_prefilt_fasta_gz  = ch_assemblies_prefilt_fasta_gz.mix(ch_penguin_co_fasta_gz)
        }

        // REMOVE EMPTY FASTA FILES FROM CHANNEL
        ch_assemblies_fasta_gz = ch_assemblies_prefilt_fasta_gz
            .filter { meta, fasta ->
                try {
                    fasta.countFasta(limit: 5) > 1
                } catch (EOFException) {
                    log.warn "[HoffLab/phageannotator]: ${fasta} has an EOFException, this is likely an empty gzipped file."
                }
            }

        // REMOVE EMPTY GRAPH FILES FROM CHANNEL
        ch_assembly_graphs_gz = ch_assemblies_prefilt_fasta_gz
            .filter { meta, graph ->
                try {
                    graph.countLines(limit: 5) > 1
                } catch (EOFException) {
                    log.warn "[HoffLab/phageannotator]: ${graph} has an EOFException, this is likely an empty gzipped file."
                }
            }
    } else {
        ch_assemblies_fasta_gz  = []
        ch_assembly_graphs_gz   = []
        ch_assembly_log         = []
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
        Resolve assembly graph
    ------------------------------------------------------------------------------*/
    if (params.run_phables) {
        // create channel from params.phables_db
        if (!params.phables_db){
            ch_phables_db = null
        } else {
            ch_phables_db = Channel.value(
                file(params.phables_db, checkIfExists:true)
            )
        }

        // filter to only megahit and metaspades assemblies
        ch_assembly_graphs_phables = ch_assembly_graphs_gz
            .filter { meta, graph ->
                meta.assembler.contains('metaspades_single') || meta.assembler.contains('megahit_single')
            }

        // identify megahit graphs and convert to gfa format
        ch_branched_graphs = ch_assembly_graphs_phables
            .branch { meta, graph ->
                megahit:    meta.assembler.contains('megahit')
                metaspades: true
            }

        //
        // MODULE: Convert FastG to GFA file
        //
        FASTG2GFA(
            ch_branched_graphs.megahit
        ).gfa
        ch_megahit_gfa_gz = FASTG2GFA.out.gfa

        // combine metaspades gfa files with megahit gfa
        ch_gfa_gz = ch_branched_graphs.metaspades.mix(ch_megahit_gfa_gz)

        //
        // SUBWORKFLOW: Download and run phables
        //
        FASTQGFA_VIRUSEXTENSION_PHABLES(
            ch_bt2_fastq_gz,
            ch_gfa_gz,
            ch_phables_db
        )
        ch_phables_prefilt_fasta_gz = FASTQGFA_VIRUSEXTENSION_PHABLES.out.phables_fasta
        ch_versions                 = ch_versions.mix(FASTQGFA_VIRUSEXTENSION_PHABLES.out.versions)

        // REMOVE EMPTY FASTAS FROM CHANNEL
        ch_phables_fasta_gz = ch_phables_prefilt_fasta_gz
            .filter { meta, fasta ->
                try {
                    fasta.countFasta(limit: 5) > 1
                } catch (EOFException) {
                    log.warn "[HoffLab/phageannotator]: ${fasta} has an EOFException, this is likely an empty gzipped file."
                }
            }
    } else {
        ch_phables_fasta_gz = []
    }


    /*----------------------------------------------------------------------------
        Run quick quality analysis of assemblies
    ------------------------------------------------------------------------------*/
    if (params.run_seqkit_stats) {
        //
        // MODULE: Perform quick stats on assemblies
        //
        SEQKIT_STATS(
            ch_phables_fasta_gz
        )
        ch_seqkit_stats_tsv = SEQKIT_STATS.out.stats
        ch_versions         = ch_versions.mix(SEQKIT_STATS.out.versions)
    } else {
        ch_seqkit_stats_tsv = []
    }


    /*----------------------------------------------------------------------------
        Remove low-length assemblies
    ------------------------------------------------------------------------------*/
    if (params.run_seqkit_seq) {
        // combine input assemblies with generated assemblies
        ch_combined_fasta_gz = ch_raw_fasta_gz
            .mix(ch_assemblies_fasta_gz)
            .mix(ch_phables_fasta_gz)
        //
        // MODULE: Filter assemblies by length
        //
        SEQKIT_SEQ(
            ch_combined_fasta_gz
        )
        ch_seqkit_seq_prefilt_fasta_gz  = SEQKIT_SEQ.out.fastx
        ch_versions                     = ch_versions.mix(SEQKIT_SEQ.out.versions)

        // REMOVE EMPTY FASTA FILES FROM CHANNEL
        ch_seqkit_seq_prefilt_fasta_gz
            .filter { meta, fasta ->
                try {
                    fasta.countFasta(limit: 5) > 1
                } catch (EOFException) {
                    log.warn "[HoffLab/phageannotator]: ${fasta} has an EOFException, this is likely an empty gzipped file."
                }
            }
            .set { ch_seqkit_seq_fasta_gz }

        if (!params.run_cobra) {
            // identify workDirs to clean
            ch_pre_seqkit_seq_workdirs  = getWorkDirs (
                ch_assemblies_fasta_gz,
                ch_seqkit_seq_fasta_gz,
                []
            )
            ch_workdirs_to_clean        = ch_workdirs_to_clean.mix(ch_pre_seqkit_seq_workdirs)
        }
    } else {
        ch_seqkit_seq_fasta_gz  = ch_assemblies_fasta_gz
    }


    emit:
    assemblies_fasta_gz = ch_seqkit_seq_fasta_gz    // channel: [ [ meta.id, meta.group, meta.assembler ], assembly.fasta.gz ]
    multiqc_files       = ch_multiqc_files          // channel: /path.to/multiqc_files
    versions            = ch_versions               // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
