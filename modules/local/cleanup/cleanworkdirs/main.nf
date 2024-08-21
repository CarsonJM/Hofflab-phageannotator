process CLEANWORKDIRS {
    tag "${meta.id}|${module_being_cleaned}"
    label "process_single"

    input:
    tuple val(meta), val(directory), val(module_being_cleaned)

    output:
    tuple val(meta), path("*.log")  , emit: cleaned_workdir

    script:
    def files_to_clean = ""
        nextflow_files = [".command.begin", ".command.err", ".command.log", ".command.out", ".command.run", ".command.sh", ".command.trace", ".exitcode", "versions.yml"]
        file(directory, checkIfExists:true).eachFileRecurse { file ->
            files_to_clean = nextflow_files.any { file.name =~ it } ? files_to_clean : files_to_clean + " " + file
        }
    """
    for file in${files_to_clean}; do
        if ! [ -L \$file ]; then
            clean_work_files.sh "\$file" "null" >> ${meta.id}.${module_being_cleaned}.log 2>&1
        fi
    done
    """

    stub:
    def files_to_clean = ""
        nextflow_files = [".command.begin", ".command.err", ".command.log", ".command.out", ".command.run", ".command.sh", ".command.trace", ".exitcode", "versions.yml"]
        file( directory, checkIfExists:true ).eachFileRecurse { file ->
            files_to_clean = nextflow_files.any { file.name =~ it } ? files_to_clean : files_to_clean + " " + file
        }
    """
    touch ${meta.id}.${module_being_cleaned}.log
    """
}
