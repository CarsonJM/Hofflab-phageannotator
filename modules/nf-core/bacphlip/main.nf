process BACPHLIP {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-e16bfb0f667f2f3c236b32087aaf8c76a0cd2864:c64689d7d5c51670ff5841ec4af982edbe7aa406-0':
        'biocontainers/mulled-v2-e16bfb0f667f2f3c236b32087aaf8c76a0cd2864:c64689d7d5c51670ff5841ec4af982edbe7aa406-0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.bacphlip")         , emit: bacphlip_results
    tuple val(meta), path("*.hmmsearch.tsv")    , emit: hmmsearch_results
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def VERSION = '0.9.6' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    prefix = task.ext.prefix ?: "${meta.id}"
    fasta_name  = fasta.getName().replace(".gz", "")
    """
    if [ \$(zcat ${fasta} | grep ">" | wc -l) -gt 0 ]; then
        gunzip -c ${fasta} > ${fasta_name}
        bacphlip \\
            -i ${fasta_name} \\
            ${args}
    else
        touch ${fasta_name}.bacphlip
        touch ${fasta_name}.hmmsearch.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bacphlip: $VERSION
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION = '0.9.6' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch ${fasta}.bacphlip
    touch ${fasta}.hmmsearch.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bacphlip: $VERSION
    END_VERSIONS
    """
}
