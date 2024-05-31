//
// Assess virus quality with Checkv
//
include { CHECKV_DOWNLOADDATABASE   } from '../../../modules/nf-core/checkv/downloaddatabase/main'
include { CHECKV_ENDTOEND           } from '../../../modules/nf-core/checkv/endtoend/main'


workflow FASTA_VIRUSQUALITY_CHECKV {
    take:
    virus_fasta_gz  // [ [ meta ], fasta.gz ]   , assemblies/genomes (mandatory)
    checkv_db       // [ checkv_db ]            , CheckV database directory (optional)

    main:
    ch_versions = Channel.empty()

    // if checkv DB exists, skip CHECKV_DOWNLOADDATABASE
    if ( checkv_db ){
        ch_checkv_db = checkv_db
    } else {
        //
        // MODULE: download standard CheckV database
        //
        ch_checkv_db    = CHECKV_DOWNLOADDATABASE( ).checkv_db
        ch_versions     = ch_versions.mix( CHECKV_DOWNLOADDATABASE.out.versions )
    }

    //
    // MODULE: Run CheckV end to end
    //
    CHECKV_ENDTOEND ( virus_fasta_gz, ch_checkv_db.collect() )
    ch_versions             = ch_versions.mix ( CHECKV_ENDTOEND.out.versions )

    emit:
    quality_summary_tsv = CHECKV_ENDTOEND.out.quality_summary   // [ [ meta ], quality_summary.tsv ]    , TSV file containing quality data
    proteins_faa_gz     = CHECKV_ENDTOEND.out.proteins          // [ [ meta ], proteins.faa.gz ]        , FAA file containing protein sequences
    viruses_fna_gz      = CHECKV_ENDTOEND.out.viruses           // [ [ meta ], viruses.fna.gz ]         , FASTA file containing viral sequences
    versions            = ch_versions                           // [ versions.yml ]
}
