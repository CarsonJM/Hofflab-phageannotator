process CLEANWORKDIRS {
    tag "${meta.id}"
    label "process_single"

    input:
    tuple val(meta), val(directory), val(ignore)

    output:
    val(1), emit: IS_CLEAN

    script:
    def files_to_clean = ""
        nextflow_files = [".command.begin", ".command.err", ".command.log", ".command.out", ".command.run", ".command.sh", ".command.trace", ".exitcode", "versions.yml"]
        files_to_keep = ignore + nextflow_files
        file( directory, checkIfExists:true ).eachFileRecurse { file ->
            files_to_clean = files_to_keep.any { file.name =~ it } ? files_to_clean : files_to_clean + " " + file
        }
    """
    for file in${files_to_clean}; do
        if ! [ -L \$file ]; then
            clean_work_files.sh "\$file" "null"
        fi
    done
    """
}
