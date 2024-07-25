process PLASS_PENGUIN {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plass:5.cf8933--pl5321h6a68c12_0':
        'biocontainers/plass:5.cf8933--pl5321h6a68c12_0' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("${prefix}.contigs.fasta.gz") , emit: contigs
    tuple val(meta), path("${prefix}.log")              , emit: log
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    penguin \\
        guided_nuclassemble \\
        ${fastq} \\
        ${prefix}.contigs.fasta \\
        tmp \\
        --threads ${task.cpus} \\
        ${args} >> ${prefix}.log 2>&1

        gzip ${prefix}.contigs.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        penguin: \$(echo \$(penguin 2>&1) | sed -n 's/^.*PenguiN Version: //; 3p' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}.contigs.fasta.gz
    touch ${prefix}.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        penguin: \$(echo \$(penguin 2>&1) | sed -n 's/^.*PenguiN Version: //; 3p' )
    END_VERSIONS
    """
}
