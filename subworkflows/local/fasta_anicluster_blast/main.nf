//
// Compare sequences by performing an all-v-all BLAST
//
include { BLAST_MAKEBLASTDB                 } from '../../../modules/nf-core/blast/makeblastdb/main'
include { BLAST_BLASTN                      } from '../../../modules/nf-core/blast/blastn/main'
include { CHECKV_ANICALC                    } from '../../../modules/local/checkv/anicalc/main'
include { CHECKV_ANICLUST                   } from '../../../modules/local/checkv/aniclust/main'
include { CHECKV_ANICLUST as CHECKV_DEREP   } from '../../../modules/local/checkv/aniclust/main'
include { SEQKIT_GREP                       } from '../../../modules/nf-core/seqkit/grep/main'
include { SEQKIT_GREP as SEQKIT_GREP_DEREP  } from '../../../modules/nf-core/seqkit/grep/main'

workflow FASTA_ANICLUSTER_BLAST {

    take:
    fasta_gz    // [ [ meta ], fasta.gz ]   , assemblies/genomes (mandatory)
    run_derep   // Boolean                  , run dereplication (optional)
    run_cluster // Boolean                  , run clustering (optional)

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Make BLASTN database
    //
    ch_blast_db = BLAST_MAKEBLASTDB ( fasta_gz ).db
    ch_versions = ch_versions.mix( BLAST_MAKEBLASTDB.out.versions )

    // join blast and fasta for input
    ch_blastn_input = fasta_gz
        .join(ch_blast_db)
        .multiMap { it ->
            fasta:  [ it[0], it[1] ]
            db:     [ it[0], it[2] ]
        }

    //
    // MODULE: Perform BLASTN all-v-all alignment
    //
    ch_blast_txt    = BLAST_BLASTN ( ch_blastn_input.fasta , ch_blastn_input.db ).txt
    ch_versions     = ch_versions = ch_versions.mix( BLAST_MAKEBLASTDB.out.versions )

    //
    // MODULE: Calculate average nucleotide identity (ANI) and alignment fraction (AF) based on BLAST
    //
    ch_blast_ani_tsv    = CHECKV_ANICALC ( ch_blast_txt ).blast_ani
    ch_versions         = ch_versions.mix( CHECKV_ANICALC.out.versions )

    // create input for CHECKV_ANICALC
    ch_aniclust_input = fasta_gz
        .join( ch_blast_ani_tsv )
        .multiMap { it ->
            fasta: [ it[0], it[1] ]
            tsv:   [ it[0], it[2] ]
        }

    if ( run_cluster ) {
        //
        // MODULE: Cluster virus sequences based on ANI and AF
        //
        ch_clusters_tsv = CHECKV_ANICLUST (
            ch_aniclust_input.fasta,
            ch_aniclust_input.tsv
        ).clusters
        ch_reps_tsv     = CHECKV_ANICLUST.out.reps
        ch_versions     = ch_versions.mix( CHECKV_ANICLUST.out.versions )

        // create input for extracting cluster representatives
        ch_seqkit_input = fasta_gz
            .join( ch_reps_tsv )
            .multiMap { it ->
                fasta: [ it[0], it[1] ]
                tsv:   [ it[2] ]
            }

        //
        // MODULE: Extract cluster representatives
        //
        ch_votu_reps_fasta_gz   = SEQKIT_GREP ( ch_seqkit_input.fasta, ch_seqkit_input.tsv ).filter
        ch_versions             = ch_versions.mix( SEQKIT_GREP.out.versions )
    } else {
        ch_votu_reps_fasta_gz = []
    }

    if ( run_derep ) {
        //
        // MODULE: Dereplicate virus sequences based on ANI and AF
        //
        ch_derep_clusters_tsv   = CHECKV_DEREP (
            ch_aniclust_input.fasta,
            ch_aniclust_input.tsv
        ).clusters
        ch_derep_reps_tsv       = CHECKV_DEREP.out.reps
        ch_versions             = ch_versions.mix( CHECKV_DEREP.out.versions )

        // create input for extracting cluster representatives
        ch_seqkit_derep_input = fasta_gz
            .join( ch_derep_reps_tsv )
            .multiMap { it ->
                fasta: [ it[0], it[1] ]
                tsv:   [ it[2] ]
            }

        //
        // MODULE: Extract dereplicated representatives
        //
        ch_derep_reps_fasta_gz  = SEQKIT_GREP_DEREP ( ch_seqkit_input.fasta, ch_seqkit_input.tsv ).filter
        ch_versions             = ch_versions.mix( SEQKIT_GREP_DEREP.out.versions )
    } else {
        ch_derep_reps_fasta_gz = []
    }


    emit:
    votu_reps       = ch_votu_reps_fasta_gz     // [ [ meta ], fasta.gz ]   , FASTA file containing votu reps
    votu_reps_tsv   = ch_reps_tsv               // [ [ meta ], tsv ]        , TSV file containing votu reps
    derep_reps      = ch_derep_reps_fasta_gz    // [ [ meta ], fasta.gz ]   , FASTA file containing dereplicated reps
    versions        = ch_versions               // [ versions.yml ]
}
