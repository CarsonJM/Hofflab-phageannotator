process CHECKV_ENDTOEND {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/checkv:1.0.1--pyhdfd78af_0':
        'biocontainers/checkv:1.0.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path db

    output:
    tuple val(meta), path ("${prefix}_quality_summary.tsv") , emit: quality_summary
    tuple val(meta), path ("${prefix}_contamination.tsv")   , emit: contamination
    tuple val(meta), path ("${prefix}_completeness.tsv")    , emit: completeness
    tuple val(meta), path ("${prefix}_viruses.fna.gz")      , emit: viruses
    tuple val(meta), path ("${prefix}_proteins.faa.gz")     , emit: proteins
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    checkv \\
        end_to_end \\
        ${args} \\
        -t ${task.cpus} \\
        -d ${db} \\
        ${fasta} \\
        ${prefix}

    mv ${prefix}/quality_summary.tsv ${prefix}_quality_summary.tsv
    mv ${prefix}/contamination.tsv ${prefix}_contamination.tsv
    mv ${prefix}/completeness.tsv ${prefix}_completeness.tsv
    gzip ${prefix}/*.fna
    cat ${prefix}/*.fna.gz > ${prefix}_viruses.fna.gz
    mv ${prefix}/tmp/proteins.faa ${prefix}_proteins.faa
    gzip ${prefix}_proteins.faa
    rm -rf ./${prefix}/*

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkv: \$(checkv -h 2>&1 | sed -n "s/Prophage Tracer V//; s/\s.*//;  2p")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}_quality_summary.tsv
    touch ${prefix}_contamination.tsv
    touch ${prefix}_completeness.tsv
    echo "" | gzip > ${prefix}_viruses.fna.gz
    echo "" | gzip > ${prefix}_proteins.faa.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkv: \$(checkv -h 2>&1  | sed -n 's/^.*CheckV v//; s/: assessing.*//; 1p')
    END_VERSIONS
    """
}
