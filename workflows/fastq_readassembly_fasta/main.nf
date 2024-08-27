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
// MODULES: Local modules
//
// include { FASTG2GFA                             } from '../../modules/local/fastg2gfa/main'
include { PLASS_PENGUIN as PENGUIN_SINGLE       } from '../../modules/local/plass/penguin/main'
include { PLASS_PENGUIN as PENGUIN_COASSEMBLY   } from '../../modules/local/plass/penguin/main'

//
// SUBWORKFLOWS: Consisting of a mix of local and nf-core/modules
//

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

workflow FASTQ_READASSEMBLY_FASTA {

    take:
    preprocessed_fastq_gz    // channel: [ [ meta ], [ reads_1.fastq.gz, reads_2.fastq.gz ] ]

    main:
    ch_versions             = Channel.empty()
    ch_multiqc_files        = Channel.empty()
    ch_workdirs_to_clean    = Channel.empty()


    /*----------------------------------------------------------------------------
        Assemble reads into contigs/scaffolds
    ------------------------------------------------------------------------------*/
    ch_assemblies_prefilt_fasta_gz  = Channel.empty()
    ch_assembly_graph_gz            = Channel.empty()


    if (params.run_metaspades_single) {
        // prepare reads for metaspades input
        ch_metaspades_single_input = preprocessed_fastq_gz
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
        // select scaffolds or contigs depending on user selection
        if (params.metaspades_use_scaffolds) {
            ch_metaspades_single_fasta_gz   = METASPADES_SINGLE.out.scaffolds
                .join(METASPADES_SINGLE.out.min_kmer)
                .join(METASPADES_SINGLE.out.max_kmer)
                .map { meta, fasta, min_kmer, max_kmer ->
                    [ meta + [ mink: min_kmer, maxk: max_kmer ], fasta ]
                }
        } else {
            ch_metaspades_single_fasta_gz   = METASPADES_SINGLE.out.contigs
                .join(METASPADES_SINGLE.out.min_kmer)
                .join(METASPADES_SINGLE.out.max_kmer)
                .map { meta, fasta, min_kmer, max_kmer ->
                    [ meta + [ mink: min_kmer, maxk: max_kmer ], fasta ]
                }
        }
        ch_assemblies_prefilt_fasta_gz      = ch_assemblies_prefilt_fasta_gz.mix(ch_metaspades_single_fasta_gz)
        ch_assembly_graph_gz                = ch_assembly_graph_gz.mix(METASPADES_SINGLE.out.gfa)
        ch_versions                         = ch_versions.mix(METASPADES_SINGLE.out.versions)
    }

    if (params.run_megahit_single) {
        //
        // MODULE: Assemble reads individually with MEGAHIT
        //
        MEGAHIT_SINGLE(
            preprocessed_fastq_gz.map { meta, fasta -> [ meta + [ assembler: 'megahit_single' ], fasta ] }
        )
        ch_megahit_single_fasta_gz      = MEGAHIT_SINGLE.out.contigs
            .join(MEGAHIT_SINGLE.out.min_kmer)
            .join(MEGAHIT_SINGLE.out.max_kmer)
            .map { meta, fasta, min_kmer, max_kmer ->
                [ meta + [ mink: min_kmer, maxk: max_kmer ], fasta ]
            }
        ch_assemblies_prefilt_fasta_gz  = ch_assemblies_prefilt_fasta_gz.mix(ch_megahit_single_fasta_gz)
        ch_assembly_graph_gz            = ch_assembly_graph_gz.mix(MEGAHIT_SINGLE.out.graph)
        ch_versions                     = ch_versions.mix(MEGAHIT_SINGLE.out.versions)
    }

    if (params.run_penguin_single) {
        //
        // MODULE: Assemble reads individually with PenguiN
        //
        PENGUIN_SINGLE(
            preprocessed_fastq_gz.map { meta, fasta -> [ meta + [ assembler: 'penguin_single' ], fasta ] }
        )
        ch_penguin_single_fasta_gz      = PENGUIN_SINGLE.out.contigs
        ch_assemblies_prefilt_fasta_gz  = ch_assemblies_prefilt_fasta_gz.mix(ch_penguin_single_fasta_gz)
        ch_versions                     = ch_versions.mix(PENGUIN_SINGLE.out.versions)
    }

    if ( params.run_metaspades_coassembly || params.run_megahit_coassembly || params.run_penguin_coassembly ) {
        // group and set group as new id
        ch_cat_coassembly_fastq_gz = preprocessed_fastq_gz
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

        //
        // MODULE: Combine reads within groups for coassembly
        //
        CAT_COASSEMBLY(
            ch_cat_coassembly_fastq_gz.coassembly
        )
        ch_coassembly_fastq_gz  = CAT_COASSEMBLY.out.reads
        ch_versions             = ch_versions.mix(CAT_COASSEMBLY.out.versions)
    }

    if ( params.run_metaspades_coassembly ) {
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
                .join(METASPADES_COASSEMBLY.out.min_kmer)
                .join(METASPADES_COASSEMBLY.out.max_kmer)
                .map { meta, fasta, min_kmer, max_kmer ->
                    [ meta + [ mink: min_kmer, maxk: max_kmer ], fasta ]
                }
        } else {
            ch_metaspades_co_fasta_gz   = METASPADES_COASSEMBLY.out.contigs
                .join(METASPADES_COASSEMBLY.out.min_kmer)
                .join(METASPADES_COASSEMBLY.out.max_kmer)
                .map { meta, fasta, min_kmer, max_kmer ->
                    [ meta + [ mink: min_kmer, maxk: max_kmer ], fasta ]
                }
        }
        ch_assemblies_prefilt_fasta_gz  = ch_assemblies_prefilt_fasta_gz.mix(ch_metaspades_co_fasta_gz)
        ch_assembly_graph_gz            = ch_assembly_graph_gz.mix(METASPADES_COASSEMBLY.out.gfa)
        ch_versions                     = ch_versions.mix(METASPADES_COASSEMBLY.out.versions)
    }

    if (params.run_megahit_coassembly) {
        //
        // MODULE: Co-assemble reads with MEGAHIT
        //
        MEGAHIT_COASSEMBLY(
            ch_coassembly_fastq_gz.map { meta, fastq -> [ meta + [ assembler: 'megahit_coassembly' ], fastq ] }
        )
        ch_megahit_co_fasta_gz          = MEGAHIT_COASSEMBLY.out.contigs
            .join(MEGAHIT_COASSEMBLY.out.min_kmer)
            .join(MEGAHIT_COASSEMBLY.out.max_kmer)
            .map { meta, fasta, min_kmer, max_kmer ->
                [ meta + [ mink: min_kmer, maxk: max_kmer ], fasta ]
            }
        ch_versions                     = ch_versions.mix(MEGAHIT_COASSEMBLY.out.versions)
        ch_assemblies_prefilt_fasta_gz  = ch_assemblies_prefilt_fasta_gz.mix(ch_megahit_co_fasta_gz)
        ch_assembly_graph_gz            = ch_assembly_graph_gz.mix(MEGAHIT_COASSEMBLY.out.graph)
    }

    if (params.run_penguin_coassembly) {
        //
        // MODULE: Co-assemble reads with PenguiN
        //
        PENGUIN_COASSEMBLY(
            ch_cat_coassembly_fastq_gz.coassembly.map { meta, fastq -> [ meta + [ assembler: 'penguin_coassembly' ], fastq ] }
        )
        ch_penguin_co_fasta_gz          = PENGUIN_COASSEMBLY.out.contigs
        ch_versions                     = ch_versions.mix(PENGUIN_COASSEMBLY.out.versions)
        ch_assemblies_prefilt_fasta_gz  = ch_assemblies_prefilt_fasta_gz.mix(ch_penguin_co_fasta_gz)
    }

    // REMOVE EMPTY FASTA FILES FROM CHANNEL
    ch_assemblies_fasta_gz = rmEmptyFastAs(ch_assemblies_prefilt_fasta_gz)


    // /*----------------------------------------------------------------------------
    //     Resolve assembly graph
    // ------------------------------------------------------------------------------*/
    // if (params.run_phables) {
    //     // create channel from params.phables_db
    //     if (!params.phables_db){
    //         ch_phables_db = null
    //     } else {
    //         ch_phables_db = Channel.value(
    //             file(params.phables_db, checkIfExists:true)
    //         )
    //     }

    //     // identify megahit graphs and convert to gfa format
    //     ch_branched_graphs = ch_assembly_graph_gz
    //         .branch { meta, graph ->
    //             megahit:    meta.assembler.contains('megahit')
    //             metaspades: meta.assembler.contains('metaspades')
    //         }

    //     //
    //     // MODULE: Convert FastG to GFA file
    //     //
    //     FASTG2GFA(
    //         ch_branched_graphs.megahit
    //     )
    //     ch_megahit_gfa_gz = FASTG2GFA.out.gfa

    //     // combine metaspades gfa files with megahit gfa
    //     ch_gfa_gz = ch_branched_graphs.metaspades.mix(ch_megahit_gfa_gz)

    //     //
    //     // SUBWORKFLOW: Download and run phables
    //     //
    //     FASTQGFA_GRAPHRESOLUTION_PHABLES(
    //         preprocessed_fastq_gz,
    //         ch_gfa_gz,
    //         ch_phables_db
    //     )
    //     ch_phables_prefilt_fasta_gz = FASTQGFA_GRAPHRESOLUTION_PHABLES.out.phables_fasta_gz
    //     ch_versions                 = ch_versions.mix(FASTQGFA_GRAPHRESOLUTION_PHABLES.out.versions)

    //     // REMOVE EMPTY FASTAS FROM CHANNEL
    //     ch_phables_fasta_gz = rmEmptyFastAs(ch_phables_prefilt_fasta_gz)
    // } else {
    //     ch_phables_fasta_gz = ch_empty_channel
    // }

    emit:
    assemblies_fasta_gz = ch_assemblies_fasta_gz    // channel: [ [ meta ], assembly.fasta.gz ]
    multiqc_files       = ch_multiqc_files          // channel: /path.to/multiqc_files
    versions            = ch_versions               // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
