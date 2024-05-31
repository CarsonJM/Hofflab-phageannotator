process PHIST {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/phist:1.0.0--py310hd03093a_0':
        'biocontainers/phist:1.0.0--py310hd03093a_0' }"

    input:
    tuple val(meta) , path(virus_fastas)
    path(bacteria_fastas)

    output:
    tuple val(meta), path("${prefix}_common_kmers.csv") , emit: common_kmers
    tuple val(meta), path("${prefix}_predictions.csv")  , emit: predictions
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    phist.py \\
        ${virus_fastas} \\
        ${bacteria_fastas} \\
        ${prefix}_common_kmers.csv \\
        ${prefix}_predictions.csv \\
        -t ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phist: \$(phist -v 2>&1 | sed '1!d; s/^.*PHIST utility //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_common_kmers.csv
    touch ${prefix}_predictions.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phist: \$(phist -v 2>&1 | sed '1!d; s/^.*PHIST utility //')
    END_VERSIONS
    """
}
