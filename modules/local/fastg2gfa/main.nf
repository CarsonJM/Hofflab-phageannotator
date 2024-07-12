process FASTG2GFA {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/phables:1.4.0--pyhdfd78af_0':
        'biocontainers/phables:1.4.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fastg)

    output:
    tuple val(meta), path("${prefix}.gfa.gz")   , emit: gfa
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    ${projectDir}/bin/gfa1/misc/fastg2gfa \\
        ${fastg} > ${prefix}.gfa

    gzip ${prefix}.gfa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phables: \$(echo \$(phables -v 2>&1) | sed -n 's/phables, version //' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.gfa.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phables: \$(echo \$(phables -v 2>&1) | sed -n 's/phables, version //' ))
    END_VERSIONS
    """
}
