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
include { COVERM_CONTIG                         } from '../../modules/local/coverm/contig/main'
include { EXTRACTVIRALASSEMBLIES                } from '../../modules/local/extractviralassemblies/main'
include { FASTG2GFA                             } from '../../modules/local/fastg2gfa/main'
include { FILTERSEQUENCES                       } from '../../modules/local/filtersequences/main'
include { LOGAN_DOWNLOAD                        } from '../../modules/local/logan/download/main'
include { NUCLEOTIDESTATS                       } from '../../modules/local/nucleotidestats/main'
include { PLASS_PENGUIN as PENGUIN_SINGLE       } from '../../modules/local/plass/penguin/main'
include { PLASS_PENGUIN as PENGUIN_COASSEMBLY   } from '../../modules/local/plass/penguin/main'
include { SEQKIT_SEQ                            } from '../../modules/local/seqkit/seq/main'
include { SEQKIT_SPLIT                          } from '../../modules/local/seqkit/split/main'
include { SEQKIT_STATS                          } from '../../modules/local/seqkit/stats/main'
include { TANTAN                                } from '../../modules/local/tantan/main'

//
// SUBWORKFLOWS: Consisting of a mix of local and nf-core/modules
//
include { CSV_SRADOWNLOAD_FETCHNGS          } from '../../subworkflows/local/csv_sradownload_fetchngs/main'
include { FASTA_ANICLUSTER_BLAST            } from '../../subworkflows/local/fasta_anicluster_blast/main'
include { FASTA_PHAGEFUNCTION_PHAROKKA      } from '../../subworkflows/local/fasta_phagefunction_pharokka/main'
include { FASTA_PHAGEHOST_IPHOP             } from '../../subworkflows/local/fasta_phagehost_iphop/main'
include { FASTA_VIRUSCLASSIFICATION_GENOMAD } from '../../subworkflows/local/fasta_virusclassification_genomad/main'
include { FASTA_VIRUSQUALITY_CHECKV         } from '../../subworkflows/local/fasta_virusquality_checkv/main'
include { FASTQ_HOSTREMOVAL_BOWTIE2         } from '../../subworkflows/local/fastq_hostremoval_bowtie2/main'
include { FASTQ_VIRUSENRICHMENT_VIROMEQC    } from '../../subworkflows/local/fastq_virusenrichment_viromeqc/main'
include { FASTQFASTA_VIRUSEXTENSION_COBRA   } from '../../subworkflows/local/fastqfasta_virusextension_cobra/main'
include { FASTQGFA_VIRUSEXTENSION_PHABLES   } from '../../subworkflows/local/fastqgfa_virusextension_phables/main'
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
include { BACPHLIP                              } from '../../modules/nf-core/bacphlip/main'
include { CAT_FASTQ as CAT_RUNMERGE             } from '../../modules/nf-core/cat/fastq/main'
include { CAT_FASTQ as CAT_COASSEMBLY           } from '../../modules/nf-core/cat/fastq/main'
include { CAT_CAT as CAT_VIRUSES                } from '../../modules/nf-core/cat/cat/main'
include { FASTP                                 } from '../../modules/nf-core/fastp/main'
include { FASTP as FASTP_PREPROCESSED           } from '../../modules/nf-core/fastp/main'
include { FASTQC as FASTQC_RAW                  } from '../../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_PREPROCESSED         } from '../../modules/nf-core/fastqc/main'
include { GENOMAD_ENDTOEND as GENOMAD_COBRA     } from '../../modules/nf-core/genomad/endtoend/main'
include { GENOMAD_ENDTOEND as GENOMAD_PHABLES   } from '../../modules/nf-core/genomad/endtoend/main'
include { MEGAHIT as MEGAHIT_SINGLE             } from '../../modules/nf-core/megahit/main'
include { MEGAHIT as MEGAHIT_COASSEMBLY         } from '../../modules/nf-core/megahit/main'
include { MULTIQC                               } from '../../modules/nf-core/multiqc/main'
include { SPADES as METASPADES_COASSEMBLY       } from '../../modules/nf-core/spades/main'
include { SPADES as METASPADES_SINGLE           } from '../../modules/nf-core/spades/main'

//
// SUBWORKFLOWS: Installed directly from nf-core/modules
//
include { paramsSummaryMultiqc      } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML    } from '../../subworkflows/nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PHAGEANNOTATOR {

    take:
    ch_fastq_gz // channel: [ [ meta.id, meta.run, meta.group ], [ reads_1.fastq.gz, reads_2.fastq.gz ] ]
    ch_fasta_gz // channel: [ [ meta.id, meta.run, meta.group ], fasta.gz ]

    main:
    ch_versions             = Channel.empty()
    ch_multiqc_files        = Channel.empty()
    ch_workdirs_to_clean    = Channel.empty()


    // /*----------------------------------------------------------------------------
    //     Read download
    // ------------------------------------------------------------------------------*/
    if ( params.sra_accessions ) {
        ch_sra_fastq    = CSV_SRADOWNLOAD_FETCHNGS( params.sra_accessions, params.sra_download_method ).fastq
        ch_versions     = ch_versions.mix( CSV_SRADOWNLOAD_FETCHNGS.out.versions )
        ch_fastq_gz     = ch_fastq_gz.mix( ch_sra_fastq )
    }

    //
    // MODULE: Run FastQC on raw reads
    //
    FASTQC_RAW (
        ch_fastq_gz
    )
    ch_multiqc_files    = ch_multiqc_files.mix( FASTQC_RAW.out.zip.collect{it[1]} )
    ch_versions         = ch_versions.mix( FASTQC_RAW.out.versions.first() )


    /*----------------------------------------------------------------------------
        Read merging
    ------------------------------------------------------------------------------*/
    if ( params.perform_run_merging ) {
        // prepare reads for concatenating within runs
        ch_reads_forcat             = ch_fastq_gz
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
        ch_merged_reads_fastq_gz    = ch_cat_reads_fastq_gz.mix( ch_reads_forcat_skipped )
        ch_versions                 = ch_versions.mix( CAT_RUNMERGE.out.versions )

        // identify workDirs to clean
        ch_pre_merge_workdirs       = getWorkDirs(
            ch_fastq_gz,
            ch_merged_reads_fastq_gz,
            []
        )
        ch_workdirs_to_clean        = ch_workdirs_to_clean.mix( ch_pre_merge_workdirs )
    } else {
        ch_merged_fastq_gz = ch_fastq_gz
    }


    /*----------------------------------------------------------------------------
        Read Preprocessing
    ------------------------------------------------------------------------------*/
    if ( params.run_fastp ) {
        //
        // MODULE: Run fastp on raw reads
        //
        ch_fastp_prefilt_fastq_gz   = FASTP(
            ch_merged_reads_fastq_gz,
            [],
            false,
            false
        ).reads
        ch_versions         = ch_versions.mix( FASTP.out.versions )
        ch_multiqc_files    = ch_multiqc_files.mix( FASTP.out.json.collect{ it[1] } )

        // remove empty fastQ files from channel
        ch_fastp_fastq_gz = ch_fastp_prefilt_fastq_gz.filter { meta, fastq ->
            if ( meta.single_end ) {
                try {
                    fastq.countFastq( limit: 10 ) > 1
                } catch (EOFException) {
                    log.warn "[HoffLab/phageannotator]: ${fastq} has an EOFException, this is likely an empty gzipped file."
                }
            } else {
                try {
                    fastq[0].countLines( limit: 10 ) > 1
                } catch (EOFException) {
                    log.warn "[HoffLab/phageannotator]: ${fastq[0]} has an EOFException, this is likely an empty gzipped file."
                }
            }
        }

        // identify workDirs to clean
        ch_pre_fastp_workdirs   = getWorkDirs(
            ch_merged_reads_fastq_gz,
            ch_fastp_fastq_gz,
            []
        )
        ch_workdirs_to_clean    = ch_workdirs_to_clean.mix( ch_pre_fastp_workdirs )
    } else {
        ch_fastp_fastq_gz = ch_merged_fastq_gz
    }


    /*----------------------------------------------------------------------------
        Host read removal
    ------------------------------------------------------------------------------*/
    if ( params.run_bowtie2_host_removal ) {
        // prepare host fasta and bowtie2 index
        if (params.bowtie2_igenomes_host) {
            ch_bowtie2_host_fasta   = Channel.value(
                [ [ id:'bowtie2_fasta' ], file( params.genomes[params.bowtie2_igenomes_host].fasta, checkIfExists: true ) ]
            )
            ch_bowtie2_host_index   = Channel.value(
                [ [ id:'bowtie2_index' ], file( params.genomes[params.bowtie2_igenomes_host].bowtie2, checkIfExists: true ) ]
            )
        } else {
            ch_bowtie2_host_fasta   = Channel.value(
                [ [ id:'bowtie2_fasta' ], file( params.bowtie2_custom_host_fasta, checkIfExists: true ) ]
            )
            ch_bowtie2_host_index   = null
        }

        //
        // SUBWORKFLOW: Remove host reads using Bowtie2
        //
        ch_bt2_prefilt_fastq_gz = FASTQ_HOSTREMOVAL_BOWTIE2(
            ch_fastp_fastq_gz,
            ch_bowtie2_host_fasta,
            ch_bowtie2_host_index
        ).fastq_gz
        ch_versions         = ch_versions.mix( FASTQ_HOSTREMOVAL_BOWTIE2.out.versions )
        ch_multiqc_files    = ch_multiqc_files.mix( FASTQ_HOSTREMOVAL_BOWTIE2.out.mqc.collect{ it[1] } )

        // remove empty fastQ files from channel
        ch_bt2_fastq_gz = ch_bt2_prefilt_fastq_gz.filter { meta, fastq ->
            if ( meta.single_end ) {
                try {
                    fastq.countFastq( limit: 10 ) > 1
                } catch (EOFException) {
                    log.warn "[HoffLab/phageannotator]: ${fastq} has an EOFException, this is likely an empty gzipped file."
                }
            } else {
                try {
                    fastq[0].countLines( limit: 10 ) > 1
                } catch (EOFException) {
                    log.warn "[HoffLab/phageannotator]: ${fastq[0]} has an EOFException, this is likely an empty gzipped file."
                }
            }
        }

        // identify workDirs to clean
        ch_pre_bt2_workdirs     = getWorkDirs(
            ch_fastp_fastq_gz,
            ch_bt2_fastq_gz,
            []
        )
        ch_workdirs_to_clean    = ch_workdirs_to_clean.mix( ch_pre_bt2_workdirs )
    } else {
        ch_bt2_fastq_gz  = ch_fastp_fastq_gz
    }


    /*----------------------------------------------------------------------------
        Preprocessing analysis
    ------------------------------------------------------------------------------*/
    if ( params.perform_run_merging || params.run_fastp || params.run_bowtie2_host_removal ) {
        //
        // MODULE: Run FastQC on preprocessed reads
        //
        FASTQC_PREPROCESSED (
            ch_bt2_fastq_gz
        )
        ch_multiqc_files    = ch_multiqc_files.mix( FASTQC_PREPROCESSED.out.zip.collect{it[1]} )
        ch_versions         = ch_versions.mix( FASTQC_PREPROCESSED.out.versions.first() )

        //
        // MODULE: Run fastp on preprocessed reads
        //
        FASTP_PREPROCESSED (
            ch_bt2_fastq_gz,
            [],
            false,
            false
        )
        ch_multiqc_files    = ch_multiqc_files.mix( FASTP_PREPROCESSED.out.json.collect{it[1]} )
        ch_versions         = ch_versions.mix( FASTQC_PREPROCESSED.out.versions.first() )
    }



    /*----------------------------------------------------------------------------
        Estimate virus enrichment
    ------------------------------------------------------------------------------*/
    if ( params.run_viromeqc ) {
        // create channel from params.viromeqc_db
        if ( !params.viromeqc_db ) {
            ch_viromeqc_db  = null
        } else {
            ch_viromeqc_db  = Channel.value(
                file( params.viromeqc_db, checkIfExists:true )
            )
        }

        //
        // SUBWORKFLOW: Estimate virus enrichment with ViromeQC
        //
        ch_vqc_enrich_tsv   = FASTQ_VIRUSENRICHMENT_VIROMEQC( ch_bt2_fastq_gz, ch_viromeqc_db ).enrichment_tsv
        ch_versions         = ch_versions.mix( FASTQ_VIRUSENRICHMENT_VIROMEQC.out.versions )
    } else {
        ch_vqc_enrich_tsv   = []
    }


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

        if ( params.run_metaspades_single ) {
            // prepare reads for metaspades input
            ch_metaspades_single_input  = ch_bt2_fastq_gz
                .map { meta, fastq ->
                    [ meta, fastq, [], [] ]
                }
                .branch { meta, fastq, extra1, extra2  ->
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
            if ( params.metaspades_use_scaffolds ) {
                ch_metaspades_single_fasta_gz   = METASPADES_SINGLE.out.scaffolds.map { meta, fasta ->
                    [ meta + [ assembler: 'metaspades_single', coassembly: false ], fasta ]
                }
            } else {
                ch_metaspades_single_fasta_gz   = METASPADES_SINGLE.out.contigs.map { meta, fasta ->
                    [ meta + [ assembler: 'metaspades_single', coassembly: false ], fasta ]
                }
            }
            ch_assemblies_prefilt_fasta_gz  = ch_assemblies_prefilt_fasta_gz.mix( ch_metaspades_single_fasta_gz )
            ch_assembly_graph_prefilt_gz    = ch_assembly_graph_gz.mix(
                METASPADES_SINGLE.out.gfa.map { meta, graph ->
                    [ meta + [ assembler: 'metaspades_single', coassembly: false ], graph ]
                }
            )
            ch_versions                     = ch_versions.mix( METASPADES_SINGLE.out.versions )
        }

        if ( params.run_megahit_single ) {
            //
            // MODULE: Assemble reads individually with MEGAHIT
            //
            ch_megahit_single_fasta_gz  = MEGAHIT_SINGLE( ch_bt2_fastq_gz ).contigs.map { meta, fasta ->
                [ meta + [ assembler: 'megahit_single', coassembly: false ], fasta ]
            }

            ch_assemblies_prefilt_fasta_gz  = ch_assemblies_prefilt_fasta_gz.mix( ch_megahit_single_fasta_gz )
            ch_assembly_graph_prefilt_gz    = ch_assembly_graph_gz.mix(
                MEGAHIT_SINGLE.out.graph.map { meta, graph ->
                    [ meta + [ assembler: 'megahit_single', coassembly: false ], graph ]
                }
            )
            ch_versions                     = ch_versions.mix( MEGAHIT_SINGLE.out.versions )
        }

        if ( params.run_penguin_single ) {
            //
            // MODULE: Assemble reads individually with PenguiN
            //
            ch_penguin_single_fasta_gz  = PENGUIN_SINGLE( ch_bt2_fastq_gz ).contigs.map { meta, fasta ->
                [ meta + [ assembler: 'penguin_single', coassembly: false ], fasta ]
            }

            ch_assemblies_prefilt_fasta_gz  = ch_assemblies_prefilt_fasta_gz.mix( ch_penguin_single_fasta_gz )
            ch_versions                     = ch_versions.mix( PENGUIN_SINGLE.out.versions )
        }

        if ( params.run_metaspades_coassembly || params.run_megahit_coassembly || params.run_penguin_coassembly ) {
            // group and set group as new id
            ch_cat_coassembly_fastq_gz      = ch_bt2_fastq_gz
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
            ch_coassembly_fastq_gz  = CAT_COASSEMBLY(
                ch_cat_coassembly_fastq_gz.coassembly
            ).reads
            ch_versions             = ch_versions.mix( CAT_COASSEMBLY.out.versions )

            // prepare reads for metaspades input
            ch_metaspades_coassembly_input  = ch_coassembly_fastq_gz
                .map { meta, fastq ->
                    [ meta, fastq, [], [] ]
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
            if ( params.metaspades_use_scaffolds ) {
                ch_metaspades_co_fasta_gz   = METASPADES_COASSEMBLY.out.scaffolds.map { meta, fasta ->
                    [ meta + [ assembler: 'metaspades_coassembly', coassembly: true ], fasta ]
                }
            } else {
                ch_metaspades_co_fasta_gz   = METASPADES_COASSEMBLY.out.contigs.map { meta, fasta ->
                    [ meta + [ assembler: 'metaspades_coassembly', coassembly: true ], fasta ]
                }
            }
            ch_assemblies_prefilt_fasta_gz  = ch_assemblies_prefilt_fasta_gz.mix( ch_metaspades_co_fasta_gz )
            ch_assembly_graph_prefilt_gz    = ch_assembly_graph_gz.mix(
                METASPADES_COASSEMBLY.out.gfa.map { meta, graph ->
                    [ meta + [ assembler: 'metaspades_coassembly', coassembly: false ], graph ]
                }
            )
            ch_versions                     = ch_versions.mix( METASPADES_COASSEMBLY.out.versions )
        }

        if ( params.run_megahit_coassembly ) {
            //
            // MODULE: Co-assemble reads with MEGAHIT
            //
            ch_megahit_co_fasta_gz = MEGAHIT_COASSEMBLY( ch_cat_coassembly_fastq_gz.coassembly ).contigs.map{ meta, fasta ->
                    [ meta + [ assembler: 'megahit_coassembly', coassembly: true ], fasta ]
                }
            ch_versions             = ch_versions.mix( MEGAHIT_COASSEMBLY.out.versions )

            ch_assemblies_prefilt_fasta_gz  = ch_assemblies_prefilt_fasta_gz.mix( ch_megahit_co_fasta_gz )
            ch_assembly_graph_prefilt_gz    = ch_assembly_graph_gz.mix(
                MEGAHIT_COASSEMBLY.out.graph.map { meta, graph ->
                    [ meta + [ assembler: 'megahit_coassembly', coassembly: false ], graph ]
                }
            )
        }

        if ( params.run_penguin_coassembly ) {
            //
            // MODULE: Co-assemble reads with PenguiN
            //
            ch_penguin_co_fasta_gz  = PENGUIN_COASSEMBLY( ch_cat_coassembly_fastq_gz.coassembly ).contigs.map{ meta, fasta ->
                    [ meta + [ assembler: 'penguin_coassembly', coassembly: true ], fasta ]
                }
            ch_versions                     = ch_versions.mix( PENGUIN_COASSEMBLY.out.versions )
            ch_assemblies_prefilt_fasta_gz  = ch_assemblies_prefilt_fasta_gz.mix( ch_penguin_co_fasta_gz )
        }

        // remove empty fastA files from channel
        ch_assemblies_fasta_gz = ch_assemblies_prefilt_fasta_gz.filter { meta, fasta ->
            try {
                fasta.countFasta( limit: 5 ) > 1
            } catch (EOFException) {
                log.warn "[HoffLab/phageannotator]: ${fasta} has an EOFException, this is likely an empty gzipped file."
            }
        }

        // remove empty grahp files from channel
        ch_assembly_graphs_gz   = ch_assemblies_prefilt_fasta_gz.filter { meta, graph ->
            try {
                graph.countLines( limit: 5 ) > 1
            } catch (EOFException) {
                log.warn "[HoffLab/phageannotator]: ${graph} has an EOFException, this is likely an empty gzipped file."
            }
        }
    } else {
        ch_assemblies_fasta_gz  = ch_fasta_gz
        ch_assembly_graphs_gz   = []
    }


    /*----------------------------------------------------------------------------
        Download Logan files
    ------------------------------------------------------------------------------*/
    if ( params.logan_accessions ) {
        ch_logan_accession_files        = Channel.fromFilePairs( params.logan_accessions, size: 1 )
        .map { meta, accession_file ->
            [ [ id:"logan_acc_" + meta ], accession_file ]
        }
        ch_logan_fasta          = LOGAN_DOWNLOAD( ch_logan_accession_files ).fasta
        ch_assemblies_fasta_gz  = ch_assemblies_fasta_gz.mix( ch_logan_fasta )
        ch_versions             = ch_versions.mix( LOGAN_DOWNLOAD.out.versions )
    }


    /*----------------------------------------------------------------------------
        Run quick quality analysis of assemblies
    ------------------------------------------------------------------------------*/
    if ( params.run_seqkit_stats ) {
        //
        // MODULE: Filter assemblies by length
        //
        ch_seqkit_stats_tsv = SEQKIT_STATS( ch_assemblies_fasta_gz ).stats
        ch_versions         = ch_versions.mix( SEQKIT_STATS.out.versions )
    } else {
        ch_seqkit_stats_tsv = []
    }


    /*----------------------------------------------------------------------------
        Remove low-length assemblies
    ------------------------------------------------------------------------------*/
    if ( params.run_seqkit_seq ) {
        //
        // MODULE: Filter assemblies by length
        //
        ch_seqkit_seq_prefilt_fasta_gz  = SEQKIT_SEQ( ch_assemblies_fasta_gz ).fastx
        ch_versions             = ch_versions.mix( SEQKIT_SEQ.out.versions )

        // remove empty fastA files from channel
        ch_seqkit_seq_fasta_gz = ch_seqkit_seq_prefilt_fasta_gz.filter { meta, fasta ->
            try {
                fasta.countFasta( limit: 5 ) > 1
            } catch (EOFException) {
                log.warn "[HoffLab/phageannotator]: ${fasta} has an EOFException, this is likely an empty gzipped file."
            }
        }

        if ( !params.run_cobra ) {
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


    /*----------------------------------------------------------------------------
        De novo virus classification
    ------------------------------------------------------------------------------*/
    if ( params.run_genomad || params.run_cobra ) {
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
        ch_genomad_summary_tsv      = FASTA_VIRUSCLASSIFICATION_GENOMAD( ch_seqkit_seq_fasta_gz, ch_genomad_db ).virus_summary_tsv
        ch_genomad_prefilt_fasta_gz = FASTA_VIRUSCLASSIFICATION_GENOMAD.out.virus_fasta_gz
        ch_versions                 = ch_versions.mix( FASTA_VIRUSCLASSIFICATION_GENOMAD.out.versions )

        // remove empty fastA files from channel
        ch_genomad_fasta_gz = ch_genomad_prefilt_fasta_gz.filter { meta, fasta ->
            try {
                fasta.countFasta( limit: 5 ) > 1
            } catch (EOFException) {
                log.warn "[HoffLab/phageannotator]: ${fasta} has an EOFException, this is likely an empty gzipped file."
            }
        }

        if ( !params.run_cobra || !params.run_phables ) {
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
    if ( params.run_cobra ) {
        // filter to only single assemblies with reads available
        ch_cobra_input = ch_genomad_fasta_gz.map { meta, fasta ->
                [ meta.id, meta, fasta ]
            }
            .combine( ch_bt2_fastq_gz.map { meta, fastq -> [ meta.id, meta, fastq ] }, by:0 )
            .map {
                id, meta_fasta, fasta, meta_fastq, fastq ->
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
        ch_cobra_prefilt_fasta_gz   = FASTQFASTA_VIRUSEXTENSION_COBRA.out.extended_fasta
        ch_cobra_summary_tsv        = FASTQFASTA_VIRUSEXTENSION_COBRA.out.cobra_summary_tsv
        ch_versions                 = ch_versions.mix(FASTQFASTA_VIRUSEXTENSION_COBRA.out.versions)

        // remove empty fastA files from channel
        ch_cobra_filt_fasta_gz = ch_cobra_prefilt_fasta_gz.filter { meta, fasta ->
            try {
                fasta.countFasta( limit: 5 ) > 1
            } catch (EOFException) {
                log.warn "[HoffLab/phageannotator]: ${fasta} has an EOFException, this is likely an empty gzipped file."
            }
        }

        //
        // MODULE: Run geNomad on extended viral contigs
        //
        ch_cobra_genomad_fasta_gz = GENOMAD_COBRA(
            ch_cobra_filt_fasta_gz,
            FASTA_VIRUSCLASSIFICATION_GENOMAD.out.genomad_db.first()
        ).virus_fasta.map { meta, fasta ->
                [ meta + [ assembler: meta.assember + '_cobra' ], fasta ]
            }

        ch_genomad_summary_tsv  = ch_genomad_summary_tsv
            .mix(
                GENOMAD_COBRA.out.virus_summary.map { meta, fasta ->
                    meta.assembler = meta.assembler + "_cobra"
                    [ meta, fasta ]
                }
            )
        ch_versions             = ch_versions.mix( GENOMAD_COBRA.out.versions )

        // combine cobra extended viruses with unextended viruses
        ch_cobra_fasta_gz = ch_genomad_fasta_gz
            .mix( ch_cobra_genomad_fasta_gz )
    } else {
        ch_cobra_fasta_gz       = ch_genomad_fasta_gz
        ch_cobra_summary_tsv    = Channel.empty()
    }


    /*----------------------------------------------------------------------------
        Resolve assembly graph
    ------------------------------------------------------------------------------*/
    if ( params.run_phables ) {
        // create channel from params.phables_db
        if ( !params.phables_db ){
            ch_phables_db   = null
        } else {
            ch_phables_db   = Channel.value(
                file( params.phables_db, checkIfExists:true )
            )
        }

        // filter to only megahit and metaspades assemblies
        ch_assembly_graphs_phables = ch_assembly_graphs_gz.filter{ meta, graph ->
            meta.assembler.contains('metaspades_single') || meta.assembler.contains('megahit_single')
            }

        // identify megahit graphs and convert to gfa format
        ch_branched_graphs = ch_assembly_graphs_phables.branch { meta, graph ->
            megahit:    meta.assembler.contains('megahit')
            metaspades: true
            }

        //
        // MODULE: Convert FastG to GFA file
        //
        ch_megahit_gfa_gz = FASTG2GFA( ch_branched_graphs.megahit ).gfa

        // combine metaspades gfa files with megahit gfa
        ch_gfa_gz = ch_branched_graphs.metaspades.mix( ch_megahit_gfa_gz )

        //
        // SUBWORKFLOW: Download and run phables
        //
        ch_phables_prefilt_fasta_gz = FASTQGFA_VIRUSEXTENSION_PHABLES(
            ch_bt2_fastq_gz,
            ch_gfa_gz,
            ch_phables_db
        ).phables_fasta
        ch_versions                 = ch_versions.mix( FASTQGFA_VIRUSEXTENSION_PHABLES.out.versions )

        // remove empty fastA files from channel
        ch_phables_filt_fasta_gz = ch_phables_prefilt_fasta_gz.filter { meta, fasta ->
            try {
                fasta.countFasta( limit: 5 ) > 1
            } catch (EOFException) {
                log.warn "[HoffLab/phageannotator]: ${fasta} has an EOFException, this is likely an empty gzipped file."
            }
        }

        //
        // MODULE: Run geNomad on extended viral contigs
        //
        ch_phables_genomad_fasta_gz = GENOMAD_PHABLES(
            ch_phables_filt_fasta_gz,
            FASTA_VIRUSCLASSIFICATION_GENOMAD.out.genomad_db.first()
        ).virus_fasta.map { meta, fasta ->
                [ meta + [ assembler: meta.assember + '_phabes' ], fasta ]
            }

        ch_genomad_summary_tsv  = ch_genomad_summary_tsv
            .mix(
                GENOMAD_PHABLES.out.virus_summary.map { meta, fasta ->
                    [ meta + [ assembler: meta.assember + '_phabes' ], fasta ]
                }
            )
        ch_versions             = ch_versions.mix( GENOMAD_PHABLES.out.versions )

        // combine phables extended viruses with other viruses
        ch_phables_fasta_gz = ch_cobra_fasta_gz
            .mix( ch_phables_genomad_fasta_gz )
    } else {
        ch_phables_fasta_gz = ch_cobra_fasta_gz
    }


    /*----------------------------------------------------------------------------
        Assess virus quality and filter
    ------------------------------------------------------------------------------*/
    if ( params.run_checkv || params.run_nucleotide_stats ) {
        // create channel from params.checkv_db
        if ( !params.checkv_db ){
            ch_checkv_db    = null
        } else {
            ch_checkv_db    = Channel.value(
                file( params.checkv_db, checkIfExists:true )
            )
        }

        //
        // SUBWORKFLOW: Assess virus quality with CheckV
        //
        ch_quality_summary_tsv  = FASTA_VIRUSQUALITY_CHECKV(ch_phables_fasta_gz, ch_checkv_db).quality_summary_tsv
        ch_checkv_faa_gz        = FASTA_VIRUSQUALITY_CHECKV.out.proteins_faa_gz
        ch_checkv_fna_gz        = FASTA_VIRUSQUALITY_CHECKV.out.viruses_fna_gz
        ch_versions             = ch_versions.mix( FASTA_VIRUSQUALITY_CHECKV.out.versions )

        // identify intermediate workDirs to clean
        ch_pre_checkv_workdirs  = getWorkDirs(
            ch_cobra_fasta_gz,
            ch_quality_summary_tsv,
            []
        )
        ch_workdirs_to_clean    = ch_workdirs_to_clean.mix( ch_pre_checkv_workdirs )
    } else {
        ch_quality_summary_tsv  = []
        ch_checkv_fna_gz        = ch_cobra_fasta_gz
    }

    if ( params.run_tantan ) {
        //
        // MODULE: Identify low-complexity regions with tantan
        //
        ch_tantan_tsv   = TANTAN( ch_checkv_fna_gz ).tantan
        ch_versions     = ch_versions.mix( TANTAN.out.versions )

        // identify intermediate workDirs to clean
        ch_pre_tantan_workdirs  = getWorkDirs(
            ch_checkv_fna_gz,
            ch_tantan_tsv,
            []
        )
        ch_workdirs_to_clean    = ch_workdirs_to_clean.mix( ch_pre_tantan_workdirs )
    } else {
        ch_tantan_tsv   = []
    }

    if ( params.run_nucleotide_stats ) {
        // join fasta and proteins for nucleotide stats input
        ch_nuc_stats_input  = ch_checkv_fna_gz
            .join( ch_checkv_faa_gz )
            .multiMap { it ->
                fasta: [ it[0], it[1] ]
                proteins: [ it[0], it[2] ]
            }

        //
        // MODULE: Calculate nucleotide stats
        //
        ch_nuc_stats_tsv    = NUCLEOTIDESTATS( ch_nuc_stats_input.fasta, ch_nuc_stats_input.proteins ).nuc_stats
        ch_versions         = ch_versions.mix( NUCLEOTIDESTATS.out.versions )

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
    if ( params.run_sequence_filtering ) {
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

        if ( params.contigs_to_keep ) {
            ch_contigs_to_keep_tsv  = Channel.fromPath(
                file( params.contigs_to_keep, checkIfExists: true )
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
        ch_filt_fasta_gz = ch_filt_prefilt_fasta_gz.filter { meta, fasta ->
            try {
                fasta.countFasta( limit: 5 ) > 1
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
        ch_workdirs_to_clean    = ch_workdirs_to_clean.mix( ch_pre_filt_workdirs )
    } else {
        ch_filt_fasta_gz        = ch_checkv_fna_gz
    }


    /*----------------------------------------------------------------------------
        Dereplicate phages
    ------------------------------------------------------------------------------*/
    if ( params.run_blast_dereplication || params.run_blast_clustering || params.run_coverm ) {
        // create a channel for combining filtered viruses (sorted so output is the same for tests)
        ch_cat_viruses_input = ch_filt_fasta_gz
            .map { [ [ id:'all_samples' ], it[1] ] }
            .groupTuple( sort: 'deep' )

        //
        // MODULE: Concatenate all quality filtered viruses into one file
        //
        ch_derep_fna_gz = CAT_VIRUSES ( ch_cat_viruses_input ).file_out
        ch_versions     = ch_versions.mix( CAT_VIRUSES.out.versions )

        //
        // SUBWORKFLOW: Calculate skani all-v-all based ANI
        //
        ch_votu_reps_fasta_gz   = FASTA_ANICLUSTER_BLAST(
            ch_derep_fna_gz,
            params.run_blast_dereplication,
            params.run_blast_clustering
        ).votu_reps
        ch_versions             = ch_versions.mix( FASTA_ANICLUSTER_BLAST.out.versions )

        if ( params.run_blast_dereplication ) {
            //
            // MODULE: Split dereplicated viruses into even sized chunks
            //
            SEQKIT_SPLIT( FASTA_ANICLUSTER_BLAST.out.derep_reps )
            ch_versions             = ch_versions.mix( SEQKIT_SPLIT.out.versions )
            // flatten and create id for each split
            ch_derep_split_fasta_gz = SEQKIT_SPLIT.out.fastx
                .map { meta, files ->
                    files
                }
                .flatten()
                .map { file ->
                    [ [ id: (file.getName() =~ /(?<=all_samples_derep\.)(.*)(?=\.fa\.gz)/)[0][1] ], file ]
                }
        }
    } else {
        ch_derep_split_fasta_gz = ch_filt_fasta_gz
    }


    /*----------------------------------------------------------------------------
        Calculate vOTU abundance by read-mapping
    ------------------------------------------------------------------------------*/
    if ( params.run_coverm ) {
        //
        // MODULE: Align reads to vOTU reps
        //
        ch_alignment_results_tsv = COVERM_CONTIG ( ch_bt2_fastq_gz, ch_votu_reps_fasta_gz.first() ).alignment_results
        ch_versions = ch_versions.mix( COVERM_CONTIG.out.versions )
    }


    /*----------------------------------------------------------------------------
        Predict phage hosts
    ------------------------------------------------------------------------------*/
    if ( params.run_iphop || params.run_phist || params.run_propagate || params.run_prophage_tracer ) {
        EXTRACTVIRALASSEMBLIES
    }
    if ( params.run_iphop ) {
        // create channel from params.checkv_db
        if ( !params.iphop_db ){
            ch_iphop_db = null
        } else {
            ch_iphop_db = file( params.iphop_db, checkIfExists:true )
        }

        //
        // SUBWORKFLOW: Download database and predict phage hosts
        //
        ch_iphop_predictions_tsv    = FASTA_PHAGEHOST_IPHOP( ch_derep_split_fasta_gz, ch_iphop_db ).host_predictions_tsv
        ch_versions                 = ch_versions.mix( FASTA_PHAGEHOST_IPHOP.out.versions )
    } else {
        ch_iphop_predictions_tsv    = Channel.empty()
    }

    if ( params.run_phist ) {
        // Split phages into individual fasta files for phist
        ch_phist_fastas = PHIST_SPLIT( ch_derep_split_fasta_gz ).dir
        ch_versions     = ch_versions.mix( PHIST_SPLIT.out.versions )

        //
        // MODULE: Run PHIST to predict phage host using shared kmers
        //
        ch_phist_predictions_tsv    = PHIST( ch_phist_fastas, params.phist_bacteria_db_dir ).predictions
        ch_phist_comm_kmers_tsv     = PHIST.out.common_kmers
        ch_versions                 = ch_versions.mix( PHIST.out.versions )
    } else {
        ch_iphop_predictions_tsv    = Channel.empty()
    }

    if ( params.run_crispr_blast ) {
        // TODO: Add crispr blast option
    } else {
        ch_crispr_blast_tsv = Channel.empty()
    }


    /*----------------------------------------------------------------------------
        Phage functional annotation
    ------------------------------------------------------------------------------*/
    if ( params.run_pharokka ) {
        // create channel from params.pharokka_db
        if ( !params.pharokka_db ){
            ch_pharokka_db = null
        } else {
            ch_pharokka_db = Channel.value(
                file( params.pharokka_db, checkIfExists:true )
                )
        }

        //
        // SUBWORKFLOW: Functionally annotate phage sequences
        //
        ch_pharokka_output_tsv  = FASTA_PHAGEFUNCTION_PHAROKKA( ch_derep_split_fasta_gz, ch_pharokka_db ).pharokka_final_output_tsv
        ch_versions             = ch_versions.mix( FASTA_PHAGEFUNCTION_PHAROKKA.out.versions )
    } else {
        ch_pharokka_gbk_gz      = Channel.empty()
        ch_pharokka_output_tsv  = Channel.empty()
    }


    /*----------------------------------------------------------------------------
        Predict virus lifestyle
    ------------------------------------------------------------------------------*/
    if ( params.run_bacphlip ) {
        //
        // MODULE: Predict phage lifestyle with BACPHLIP
        //
        ch_bacphlip_lifestyle_tsv   = BACPHLIP( ch_derep_split_fasta_gz ).bacphlip_results
        ch_versions                 = ch_versions.mix( BACPHLIP.out.versions )
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
