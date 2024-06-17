process COVERM_CONTIG {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/coverm:0.7.0--h07ea13f_1':
        'biocontainers/coverm:0.7.0--h07ea13f_1' }"

    input:
    tuple val(meta), path(fastq)
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${prefix}_alignment_results.tsv")    , emit: alignment_results
    tuple val(meta), path("${prefix}.bam")                      , emit: bam
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    if [ \$(zcat ${fastq[0]} | grep "@" | wc -l) -gt 0 ]; then
        coverm \\
            contig \\
            -1 ${fastq[0]} \\
            -2 ${fastq[1]} \\
            --reference ${fasta} \\
            --output-file ${prefix}_alignment_results.tsv \\
            --bam-file-cache-directory ${prefix}_bam_files/ \\
            --threads ${task.cpus} \\
            ${args}

        mv ${prefix}_bam_files/*.bam ${prefix}.bam
    else
        touch ${prefix}_alignment_results.tsv
        touch ${prefix}.bam
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coverm: \$(echo \$(coverm --version 2>&1) | sed 's/^.*coverm //; s/Using.*\$//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}_bam_files
    touch ${prefix}_alignment_results.tsv
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        coverm: \$(echo \$(coverm --version 2>&1) | sed 's/^.*coverm //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
