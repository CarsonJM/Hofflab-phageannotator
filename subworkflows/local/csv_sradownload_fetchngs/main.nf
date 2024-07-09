//
// Download SRA metadata and FastQ files
//
// MODULES
include { ASPERACLI             } from '../../../modules/local/fetchngs/asperacli/main'
include { SRA_FASTQ_FTP         } from '../../../modules/local/fetchngs/sra_fastq_ftp/main'
include { SRA_IDS_TO_RUNINFO    } from '../../../modules/local/fetchngs/sra_ids_to_runinfo/main'
include { SRA_RUNINFO_TO_FTP    } from '../../../modules/local/fetchngs/sra_runinfo_to_ftp/main'
// SUBWORKFLOWS
include { FASTQ_DOWNLOAD_PREFETCH_FASTERQDUMP_SRATOOLS  } from '../../nf-core/fastq_download_prefetch_fasterqdump_sratools/main'


workflow CSV_SRADOWNLOAD_FETCHNGS {
    take:
    accessions          // [ accessions.tsv ]   , file containing SRA accessions (mandatory)
    sra_download_method // 'aspera_cli'         , string describing desired download method

    main:
    ch_versions = Channel.empty()

    /*----------------------------------------------------------------------------
        Load and parse accessions file
    ------------------------------------------------------------------------------*/
    //
    // Create channels from input file provided through params.input and params.assembly_input
    //
    ch_sra_accessions = Channel
        .fromSamplesheet("sra_accessions")


    /*----------------------------------------------------------------------------
        Get SRA run information
    ------------------------------------------------------------------------------*/
    //
    // MODULE: Get SRA run information for public database ids
    //
    SRA_IDS_TO_RUNINFO (
        ch_sra_accessions,
        ''
    )
    ch_versions = ch_versions.mix(SRA_IDS_TO_RUNINFO.out.versions.first())

    //
    // MODULE: Parse SRA run information, create file containing FTP links and read into workflow as [ meta, [reads] ]
    //
    SRA_RUNINFO_TO_FTP (
        SRA_IDS_TO_RUNINFO.out.tsv
    )
    ch_versions = ch_versions.mix(SRA_RUNINFO_TO_FTP.out.versions.first())

    SRA_RUNINFO_TO_FTP
        .out
        .tsv
        .splitCsv(header:true, sep:'\t')
        .map {
            meta, runinfo ->
                return [ meta + [ single_end: runinfo.single_end.toBoolean(), md5_1: runinfo.md5_1, md5_2: runinfo.md5_2 ], runinfo ]
        }
        .unique()
        .set { ch_sra_metadata }

    ch_sra_metadata
        .branch {
            meta, runinfo ->
                def download_method = 'ftp'
                // meta.fastq_aspera is a metadata string with ENA fasp links supported by Aspera
                    // For single-end: 'fasp.sra.ebi.ac.uk:/vol1/fastq/ERR116/006/ERR1160846/ERR1160846.fastq.gz'
                    // For paired-end: 'fasp.sra.ebi.ac.uk:/vol1/fastq/SRR130/020/SRR13055520/SRR13055520_1.fastq.gz;fasp.sra.ebi.ac.uk:/vol1/fastq/SRR130/020/SRR13055520/SRR13055520_2.fastq.gz'
                if ( runinfo.fastq_aspera && sra_download_method == 'aspera' ) {
                    download_method = 'aspera'
                }
                if ( ( !runinfo.fastq_aspera && !runinfo.fastq_1 ) || sra_download_method == 'sratools') {
                    download_method = 'sratools'
                }

                aspera: download_method == 'aspera'
                    return [ meta, runinfo.fastq_aspera.tokenize(';').take(2) ]
                ftp: download_method == 'ftp'
                    return [ meta, [ runinfo.fastq_1, runinfo.fastq_2 ] ]
                sratools: download_method == 'sratools'
                    return [ meta, runinfo.run_accession ]
        }
        .set { ch_sra_reads }



    /*----------------------------------------------------------------------------
        Download FastQ files
    ------------------------------------------------------------------------------*/
    if ( sra_download_method == 'ftp' ) {
        //
        // MODULE: If FTP link is provided in run information then download FastQ directly via FTP and validate with md5sums
        //
        ch_fastq_gz = SRA_FASTQ_FTP (
            ch_sra_reads.ftp
        ).fastq
        ch_versions = ch_versions.mix(SRA_FASTQ_FTP.out.versions.first())
    }

    if ( sra_download_method == 'sratools' ) {
        //
        // SUBWORKFLOW: Download sequencing reads without FTP links using sra-tools.
        //
        ch_fastq_gz = FASTQ_DOWNLOAD_PREFETCH_FASTERQDUMP_SRATOOLS (
            ch_sra_reads.sratools,
            params.dbgap_key ? file(params.dbgap_key, checkIfExists: true) : []
        ).reads
        ch_versions = ch_versions.mix(FASTQ_DOWNLOAD_PREFETCH_FASTERQDUMP_SRATOOLS.out.versions.first())
    }

    if ( sra_download_method == 'aspera' ) {
        //
        // MODULE: If Aspera link is provided in run information then download FastQ directly via Aspera CLI and validate with md5sums
        //
        ch_fastq_gz = ASPERACLI (
            ch_sra_reads.aspera,
            'era-fasp'
        ).fastq
        ch_versions = ch_versions.mix(ASPERACLI.out.versions.first())
    }

    // Isolate FASTQ channel which will be added to emit block
    ch_fastq_gz
        .map {
            meta, fastq ->
                def reads = fastq instanceof List ? fastq.flatten() : [ fastq ]
                def meta_new = meta - meta.subMap('md5_1', 'md5_2')

                fastqs = !meta.single_end ? [ reads[0], reads[1] ] : [ reads[0] ]

                return [ meta_new, fastqs ]
        }
        .set { ch_sra_metadata }

    emit:
    fastq       = ch_sra_metadata
    versions    = ch_versions.unique()
}
