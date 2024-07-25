//
// Download SRA metadata and FastQ files
//
include { fromSamplesheet       } from 'plugin/nf-validation'
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
        Get SRA run information
    ------------------------------------------------------------------------------*/
    //
    // MODULE: Get SRA run information for public database ids
    //
    SRA_IDS_TO_RUNINFO(
        accessions,
        ''
    )
    ch_versions = ch_versions.mix(SRA_IDS_TO_RUNINFO.out.versions.first())

    //
    // MODULE: Parse SRA run information, create file containing FTP links and read into workflow as [ meta, [reads] ]
    //
    SRA_RUNINFO_TO_FTP(
        SRA_IDS_TO_RUNINFO.out.tsv
    )
    ch_versions = ch_versions.mix(SRA_RUNINFO_TO_FTP.out.versions.first())

    ch_sra_metadata = SRA_RUNINFO_TO_FTP.out.tsv
        .splitCsv(header:true, sep:'\t')
        .map {
            meta ->
                def meta_clone          = meta.clone()
                    meta_clone.single_end   = meta_clone.single_end.toBoolean()
                return meta_clone
        }
        .unique()


    ch_sra_reads = ch_sra_metadata
        .branch {
            meta ->
                def download_method = 'ftp'
                // meta.fastq_aspera is a metadata string with ENA fasp links supported by Aspera
                // For single-end: 'fasp.sra.ebi.ac.uk:/vol1/fastq/ERR116/006/ERR1160846/ERR1160846.fastq.gz'
                // For paired-end: 'fasp.sra.ebi.ac.uk:/vol1/fastq/SRR130/020/SRR13055520/SRR13055520_1.fastq.gz;fasp.sra.ebi.ac.uk:/vol1/fastq/SRR130/020/SRR13055520/SRR13055520_2.fastq.gz'
                if ( meta.fastq_aspera && sra_download_method == 'aspera' ) {
                    download_method = 'aspera'
                }
                if ( ( !meta.fastq_aspera && !meta.fastq_1 ) || sra_download_method == 'sratools') {
                    download_method = 'sratools'
                }

                aspera: download_method == 'aspera'
                    return [ meta, meta.fastq_aspera.tokenize(';').take(2) ]
                ftp: download_method == 'ftp'
                    return [ meta, [ meta.fastq_1, meta.fastq_2 ] ]
                sratools: download_method == 'sratools'
                    return [ meta, meta.run_accession ]
        }


    /*----------------------------------------------------------------------------
        Download FastQ files
    ------------------------------------------------------------------------------*/
    if (sra_download_method == 'ftp') {
        //
        // MODULE: If FTP link is provided in run information then download FastQ directly via FTP and validate with md5sums
        //
        SRA_FASTQ_FTP(
            ch_sra_reads.ftp
        )
        ch_fastq_gz = SRA_FASTQ_FTP.out.fastq
        ch_versions = ch_versions.mix(SRA_FASTQ_FTP.out.versions.first())
    }

    if (sra_download_method == 'sratools') {
        //
        // SUBWORKFLOW: Download sequencing reads without FTP links using sra-tools.
        //
        FASTQ_DOWNLOAD_PREFETCH_FASTERQDUMP_SRATOOLS(
            ch_sra_reads.sratools,
            params.dbgap_key ? file(params.dbgap_key, checkIfExists: true) : []
        )
        ch_fastq_gz = FASTQ_DOWNLOAD_PREFETCH_FASTERQDUMP_SRATOOLS.out.reads
        ch_versions = ch_versions.mix(FASTQ_DOWNLOAD_PREFETCH_FASTERQDUMP_SRATOOLS.out.versions.first())
    }

    if (sra_download_method == 'aspera') {
        //
        // MODULE: If Aspera link is provided in run information then download FastQ directly via Aspera CLI and validate with md5sums
        //
        ASPERACLI(
            ch_sra_reads.aspera,
            'era-fasp'
        )
        ch_fastq_gz = ASPERACLI.out.fastq
        ch_versions = ch_versions.mix(ASPERACLI.out.versions.first())
    }

    // Isolate FASTQ channel which will be added to emit block
    ch_sra_metadata = ch_fastq_gz
        .map {
            meta, fastq ->
                def reads = fastq instanceof List ? fastq.flatten() : [ fastq ]
                def meta_new = meta - meta.subMap('md5_1', 'md5_2')

                fastqs = !meta.single_end ? [ reads[0], reads[1] ] : [ reads[0] ]

                return [ meta_new, fastqs ]
        }

    emit:
    fastq       = ch_sra_metadata
    versions    = ch_versions.unique()
}
