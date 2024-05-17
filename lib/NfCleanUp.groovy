// create a function to identify work dirs to clean
def getWorkDirs( ch_to_clean, ch_dependent, files_to_keep ) {
    // combine channel to clean and dependent channel to clean only channels that have an output
    ch_workdirs = ch_to_clean.combine ( ch_dependent, by:0 )
    // filter to retain work directory
    .map { meta, files_to_clean, dependent_files ->
        // do not clean directory if it is not a work directory
        if (( files_to_clean =~ /(^.*\/work\/[^\/]+\/[^\/]+\/).*/ )) {
            def dir_to_clean = ( files_to_clean =~ /(^.*\/work\/[^\/]+\/[^\/]+\/).*/ )[0][1]
                return [ [ id: meta.id ], dir_to_clean, files_to_keep ]
        }
    }
    // remove redundancy
    .unique()
    return ch_workdirs
}
