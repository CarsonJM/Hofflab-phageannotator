//
// Extend contig ends using COBRA
//
include { COVERM_CONTIG as COVERM_COBRA     } from '../../../modules/local/coverm/contig/main'
include { COBRAMETA                         } from '../../../modules/nf-core/cobrameta/main'

workflow FASTQFASTA_VIRUSEXTENSION_COBRA {
    take:
    fastq_gz            // [ [ meta ], [ reads_1.fastq.gz, reads_2.fastq.gz ] ] , read files (mandatory)
    fasta_gz            // [ [ meta ], fasta.gz ]                               , assemblies/scaffolds (mandatory)
    virus_summary_tsv   // [ [ meta ], virus_contigs.tsv ]                      , TSV file containing geNomad's virus summary
    cobra_assembler     // val:string {metaspades|megahit|idba}                 , string specifying which assembler was used
    mink                // val:int [ 0-999 ]                                    , integer specifying minimum kmer size used by assembler
    maxk                // val:int [ 0-999 ]                                    , integer specifying maximum kmer size used by assembler

    main:
    ch_versions = Channel.empty()

    // join fastq and fasta by meta.id
    ch_coverm_input = fastq_gz
        .join(fasta_gz)
        .multiMap { meta, fastq, fasta ->
            fastq: [ meta, fastq ]
            fasta: [ meta, fasta ]
        }

    //
    // MODULE: Calculate abundance metrics from BAM file
    //
    ch_coverage_tsv         = COVERM_COBRA(ch_coverm_input.fastq, ch_coverm_input.fasta).alignment_results
    ch_fasta_alignment_bam  = COVERM_COBRA.out.bam
    ch_versions             = ch_versions.mix(COVERM_COBRA.out.versions)

    // prepare input for cobra
    ch_cobra_input = fasta_gz
        .join(ch_coverage_tsv)
        .join(virus_summary_tsv)
        .join(ch_fasta_alignment_bam)
        .multiMap { it ->
            fasta: [ it[0], it[1] ]
            coverage: [ it[0], it[2] ]
            query: [ it[0], it[3] ]
            bam: [ it[0], it[4] ]
        }

    //
    // MODULE: Extend contigs using COBRA
    //
    ch_extended_fasta_gz = COBRAMETA (
        ch_cobra_input.fasta,
        ch_cobra_input.coverage,
        ch_cobra_input.query,
        ch_cobra_input.bam,
        cobra_assembler,
        mink,
        maxk
    ).fasta
    ch_versions         = ch_versions.mix( COBRAMETA.out.versions )

    emit:
    extended_fasta      = ch_extended_fasta_gz          // [ [ meta ], extended_contigs.fna.gz ]    , FASTA file containing extended contigs
    cobra_summary_tsv   = COBRAMETA.out.joining_summary // [ [ meta ], cobra_summary.tsv ]          , TSV file containing COBRA summary
    versions            = ch_versions.unique()          // [ versions.yml ]
}
