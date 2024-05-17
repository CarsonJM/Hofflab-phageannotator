//
// Remove host reads with bowtie2
//
include { BOWTIE2_BUILD } from '../../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN } from '../../modules/nf-core/bowtie2/build/main'

workflow FASTQ_HOSTREMOVAL_BOWTIE2 {
    take:
    fastq_gz    // channel: [ [ meta.id, meta.group ], [ reads_1.fastq.gz, reads_2.fastq.gz ] ]

    main:
    ch_versions = Channel.empty()

    // identify bowtie2 index or fasta to use for host removal
    if (params.host_genome) {
        host_fasta              = params.genomes[params.host_genome].fasta ?: false
        ch_host_fasta           = Channel
            .value(
                file("${host_fasta}", checkIfExists: true)
            )
        host_bowtie2index       = params.genomes[params.host_genome].bowtie2 ?: false
        ch_host_bowtie2index    = Channel
            .value(
                file("${host_bowtie2index}/*", checkIfExists: true)
            )
    } else if ( params.host_fasta ) {
        ch_host_fasta           = Channel
            .value(
                file("${params.host_fasta}", checkIfExists: true)
            )
    }

    //
    // MODULE: Create bowtie2 database from host FASTA
    //
    if (!ch_host_bowtie2index) {
        ch_host_bowtie2index    = BOWTIE2_BUILD (ch_host_fasta).index
        ch_versions             = ch_versions.mix(BOWTIE2_BUILD.out.versions)
    }

    //
    // MODULE: Align reads to bowtie2 index
    //
    ch_bt2_reads_fastq_gz       = BOWTIE2_ALIGN (ch_short_reads_prepped,ch_host_bowtie2index).reads
    ch_versions                 = ch_versions.mix(BOWTIE2_HOST_REMOVAL_ALIGN.out.versions)

    emit:
    bt2_reads_fastq_gz          = ch_bt2_reads_fastq_gz // channel: [ [ meta.id, meta.group ], [ reads_1.fastq.gz, reads_2.fastq.gz ] ]
    versions                    = ch_versions           // [ versions.yml ]
}
