// create a function to identify work dirs to clean
def getWorkDirs(ch_to_clean, ch_dependent) {
    // combine channel to clean and dependent channel to clean only channels that have an output
    ch_workdirs = ch_to_clean.combine ( ch_dependent, by:0 )
        // filter to retain work directory
        .map { meta, files_to_clean, dependent_files ->
            // do not clean directory if it is not a work directory
            if (( files_to_clean =~ /(^.*\/work\/[^\/]+\/[^\/]+\/).*/ )) {
                def dir_to_clean = ( files_to_clean =~ /(^.*\/work\/[^\/]+\/[^\/]+\/).*/ )[0][1]
                    return [ [ id: meta.id ], dir_to_clean ]
            }
        }
    // remove redundancy
    .unique()
    .view()
    return ch_workdirs
}

// create a function that filters channels to remove empty fasta files
def rmEmptyFastAs(ch_fastas) {
    ch_nonempty_fastas = ch_fastas
        .filter { meta, fasta ->
            try {
                fasta.countFasta(limit: 5) > 1
            } catch (java.util.zip.ZipException e) {
                log.warn "[HoffLab/phageannotator]: ${fasta} is not in GZIP format, this is likely because it was cleaned with --remove_intermediate_files"
                true
            } catch (EOFException e) {
                log.warn "[HoffLab/phageannotator]: ${fasta} has an EOFException, this likely an empty gzipped file."
            }
        }
        .view()
    return ch_nonempty_fastas
}

// create a function that filters channels to remove empty fasta files
def rmEmptyFastQs(ch_fastqs) {
    ch_nonempty_fastqs = ch_fastqs
        .filter { meta, fastq ->
                if ( meta.single_end ) {
                    try {
                        fastq.countFastq(limit: 10) > 1
                    } catch (java.util.zip.ZipException e) {
                        log.warn "[HoffLab/phageannotator]: ${fastq} is not in GZIP format, this is likely because it was cleaned with --remove_intermediate_files"
                        true
                    } catch (EOFException) {
                        log.warn "[HoffLab/phageannotator]: ${fastq} has an EOFException, this is likely an empty gzipped file."
                    }
                } else {
                    try {
                        fastq[0].countLines( limit: 10 ) > 1
                    } catch (java.util.zip.ZipException e) {
                        log.warn "[HoffLab/phageannotator]: ${fastq} is not in GZIP format, this is likely because it was cleaned with --remove_intermediate_files"
                        true
                    } catch (EOFException) {
                        log.warn "[HoffLab/phageannotator]: ${fastq[0]} has an EOFException, this is likely an empty gzipped file."
                    }
                }
            }
        .view()
    return ch_nonempty_fastqs
}
