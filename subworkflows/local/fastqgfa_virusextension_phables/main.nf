//
// Extend assembly graphs using Phables
//

include { PHABLES_INSTALL   } from '../../../modules/local/phables/install/main'
include { PHABLES_RUN       } from '../../../modules/local/phables/run/main'

workflow FASTQGFA_VIRUSEXTENSION_COBRA {
    take:
    fastq_gz            // [ [ meta ], [ reads_1.fastq.gz, reads_2.fastq.gz ] ] , read files (mandatory)
    graph_gz            // [ [ meta ], (astg|gfa).gz ]                          , assembly graph (mandatory)
    phables_db          // [ [ meta ], phables_db ]                             , Phables database (mandatory)

    main:
    ch_versions = Channel.empty()

    // if genomad_db exists, skip PHABLES_INSTALL
    if ( phables_db ){
        ch_phables_db = phables_db
    } else {
        //
        // MODULE: download geNomad database
        //
        ch_phables_db   = PHABLES_INSTALL( ).phables_db
        ch_versions     = ch_versions.mix( PHABLES_INSTALL.out.versions )
    }

    // combine fastq and gfa for phables
    ch_phables_input = fastq_gz
        .join( ch_gfa_gz )
        .multiMap { meta, fastq, gfa ->
            fastq : [ meta, fastq ]
            gfa   : [ meta, gfa ]
        }

    //
    // MODULE: Run phables to extend assemblies
    //
    ch_phables_fasta_gz = PHABLES_RUN(
        ch_phables_input.fastq,
        ch_phables_input.gfa,
        ch_phables_db.first()
    ).phables_db
    ch_versions = ch_versions.mix( PHABLES_RUN.out.versions )


    emit:
    phables_fasta       = ch_phables_fasta_gz   // [ [ meta ], extended_contigs.fna.gz ]    , FASTA file containing extended contigs
    versions            = ch_versions.unique()  // [ versions.yml ]
}
