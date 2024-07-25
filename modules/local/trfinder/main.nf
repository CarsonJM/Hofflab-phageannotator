process TRFINDER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.81':
        'biocontainers/biopython:1.81' }"

    input:
    tuple val(meta) , path(fasta)

    output:
    tuple val(meta), path("${prefix}_trfinder.tsv")         , emit: tr_stats
    tuple val(meta), path("${prefix}_trfinder.fasta.gz")    , emit: tr_fasta
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    trfinder.py \\
        --input ${fasta} \\
        --prefix ${prefix} \\
        ${args}

    gzip -f ${prefix}_trfinder.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed 's/Python //' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_trfinder.tsv
    echo "" | gzip > ${prefix}_trfinder.fasta.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed 's/Python //' )
    END_VERSIONS
    """
}
