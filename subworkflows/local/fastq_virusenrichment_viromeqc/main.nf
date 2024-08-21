//
// Estimate viral enrichment in reads
//
include { VIROMEQC_INSTALL  } from '../../../modules/local/viromeqc/install/main'
include { VIROMEQC_VIROMEQC } from '../../../modules/local/viromeqc/viromeqc/main'

workflow FASTQ_VIRUSENRICHMENT_VIROMEQC {
    take:
    fastq_gz    // [ [ meta.id ] , [ reads_1.fastq.gz, reads_2.fastq.gz ]   , reads (mandatory)
    viromeqc_db // [ viromeqc_db ]                                          , viromeqc database (optional)

    main:
    ch_versions = Channel.empty()

    // if viromeqc_db exists, skip VIROMEQC_INSTALL
    if ( viromeqc_db ){
        ch_viromeqc_db  = viromeqc_db
    } else {
        //
        // MODULE: download viromeQC database
        //
        ch_viromeqc_db  = VIROMEQC_INSTALL( ).viromeqc_index
        ch_versions     = ch_versions.mix(VIROMEQC_INSTALL.out.versions)
    }

    //
    // MODULE: Estimate viral enrichment
    //
    ch_enrichment_tsv   = VIROMEQC_VIROMEQC(fastq_gz, ch_viromeqc_db.first()).enrichment
    ch_versions         = ch_versions.mix(VIROMEQC_VIROMEQC.out.versions)

    emit:
    enrichment_tsv  = ch_enrichment_tsv // [ [ meta.id ] , enrichment.tsv ]  , TSV containing viral enrichment estimates
    versions        = ch_versions       // [ versions.yml ]
}
