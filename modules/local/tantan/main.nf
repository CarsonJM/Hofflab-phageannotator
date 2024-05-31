process TANTAN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/tantan:49--h43eeafb_0':
        'biocontainers/tantan:49--h43eeafb_0' }"

    input:
    tuple val(meta) , path(fasta)

    output:
    tuple val(meta), path("${prefix}_tantan.tsv")   , emit: tantan
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    tantan \\
        ${fasta} \\
        -f 3 \\
        ${args} > ${prefix}_tantan.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tantan: \$(tantan --version 2>&1 | sed 's/^.*tantan //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_tantan.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tantan: \$(tantan --version 2>&1 | sed 's/^.*tantan //')
    END_VERSIONS
    """
}
