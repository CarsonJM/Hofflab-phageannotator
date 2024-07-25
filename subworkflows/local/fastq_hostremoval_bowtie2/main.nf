//
// Remove host reads with bowtie2
//
include { BOWTIE2_BUILD } from '../../../modules/nf-core/bowtie2/build/main'
include { BOWTIE2_ALIGN } from '../../../modules/nf-core/bowtie2/align/main'

workflow FASTQ_HOSTREMOVAL_BOWTIE2 {
    take:
    fastq_gz        // channel: [ [ meta.id, meta.group ], [ reads_1.fastq.gz, reads_2.fastq.gz ] ]
    host_fasta_gz   // channel: [ [ id: host_fasta ], /path/to/host_fasta.fasta.gz ]
    host_bt2_index  // channel: [ [ id: host_index ], /path/to/host_index.bt2 ]

    main:
    ch_versions         = Channel.empty()
    ch_multiqc_files    = Channel.empty()

    //
    // MODULE: Create bowtie2 database from host FASTA
    //
    if (!host_bt2_index) {
        BOWTIE2_BUILD(
            host_fasta_gz
        )
        ch_host_bt2_index   = BOWTIE2_BUILD.out.index
        ch_versions         = ch_versions.mix(BOWTIE2_BUILD.out.versions)
    } else {
        ch_host_bt2_index   = host_bt2_index
    }

    //
    // MODULE: Map, generate BAM with all reads and unmapped reads in FASTQ for downstream
    //
    BOWTIE2_ALIGN(
        fastq_gz,
        ch_host_bt2_index,
        host_fasta_gz,
        true,
        true
    )
    ch_versions      = ch_versions.mix(BOWTIE2_ALIGN.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(BOWTIE2_ALIGN.out.log)

    emit:
    fastq_gz    = BOWTIE2_ALIGN.out.fastq   // channel: [ val(meta), [ reads ] ]
    versions    = ch_versions               // channel: [ versions.yml ]
    mqc         = ch_multiqc_files
}
