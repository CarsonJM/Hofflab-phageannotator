// create a function to identify work dirs to clean
def getWorkDirs(ch_to_clean, ch_dependent) {
    // combine channel to clean and dependent channel to clean only channels that have an output
    ch_double_workdir1 = Channel.empty()
    ch_double_workdir2 = Channel.empty()
    ch_branch = ch_to_clean.combine(ch_dependent, by:0)
        .branch { it ->
            double_ch: it.size() > 3
            single_ch: true
        }
    ch_double_workdir1 = ch_branch.double_ch
        .map { meta, ch_to_clean1, ch_to_clean2, ch_dep ->
            return [ [ meta ], ch_to_clean1, ch_dep ]
        }
    ch_double_workdir2 = ch_branch.double_ch
        .map { meta, ch_to_clean1, ch_to_clean2, ch_dep ->
            return [ [ meta ], ch_to_clean2, ch_dep ]
        }
    ch_workdirs = ch_double_workdir1.mix(ch_double_workdir2).mix(ch_branch.single_ch)
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
    return ch_nonempty_fastas
}

// create a function that filters channels to remove empty fasta files
def rmEmptyFastQs(ch_fastqs) {
    ch_nonempty_fastqs = ch_fastqs
        .filter { meta, fastq ->
                if ( meta.single_end ) {
                    try {
                        fastq[0].countFastq(limit: 10) > 1
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
    return ch_nonempty_fastqs
}
