//
// Classify and annotate sequences with geNomad
//
include { GENOMAD_DOWNLOAD  } from '../../../modules/nf-core/genomad/download/main'
include { GENOMAD_ENDTOEND  } from '../../../modules/nf-core/genomad/endtoend/main'

workflow FASTA_VIRUSCLASSIFICATION_GENOMAD {
    take:
    fasta_gz    // [ [ meta ], fasta.gz ]   , assemblies/genomes (mandatory)
    genomad_db  // [ genomad_db ]           , genomad database directory (optional)

    main:
    ch_versions = Channel.empty()

    // if genomad_db exists, skip GENOMAD_DOWNLOAD
    if (genomad_db){
        ch_genomad_db = genomad_db
    } else {
        //
        // MODULE: download geNomad database
        //
        ch_genomad_db   = GENOMAD_DOWNLOAD( ).genomad_db
        ch_versions     = ch_versions.mix(GENOMAD_DOWNLOAD.out.versions)
    }

    //
    // MODULE: Run geNomad end-to-end
    //
    GENOMAD_ENDTOEND(fasta_gz, ch_genomad_db)
    ch_versions = ch_versions.mix(GENOMAD_ENDTOEND.out.versions)


    emit:
    virus_summary_tsv   = GENOMAD_ENDTOEND.out.virus_summary    // [ [ meta ], summary.tsv ]    , TSV file containing genomad virus summary
    virus_fasta_gz      = GENOMAD_ENDTOEND.out.virus_fasta      // [ [ meta ], virus.fasta.gz ] , virus sequences in FASTA format
    genomad_db          = ch_genomad_db                         // [ genomad_db ]               , genomad database directory
    versions            = ch_versions                           // [ versions.yml ]
}
