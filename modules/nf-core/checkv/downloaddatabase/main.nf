process CHECKV_DOWNLOADDATABASE {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/checkv:1.0.1--pyhdfd78af_0':
        'biocontainers/checkv:1.0.1--pyhdfd78af_0' }"

    output:
    path "${prefix}/"           , emit: checkv_db
    path "${prefix}/README.txt" , emit: readme
    path "${prefix}/**.tsv"     , emit: tsv
    path "${prefix}/**.faa"     , emit: faa
    path "${prefix}/**.fna"     , emit: fna
    path "${prefix}/**.hmm"     , emit: hmm
    path "${prefix}/**.dmnd"    , emit: dmnd
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args        = task.ext.args ?: ''
    prefix          = task.ext.prefix ?: "checkv_db"

    """
    checkv download_database \\
        download \\
        ${args}

    mv download/* ${prefix}/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkv: \$(checkv -h 2>&1  | sed -n 's/^.*CheckV v//; s/: assessing.*//; 1p')
    END_VERSIONS
    """

    stub:
    def args    = task.ext.args ?: ''
    prefix      = task.ext.prefix ?: "checkv_db"

    """
    mkdir -p ${prefix}/genome_db
    touch ${prefix}/README.txt
    touch ${prefix}/genome_db/changelog.tsv
    touch ${prefix}/genome_db/checkv_error.tsv
    touch ${prefix}/genome_db/checkv_info.tsv
    touch ${prefix}/genome_db/checkv_reps.faa
    touch ${prefix}/genome_db/checkv_reps.fna
    touch ${prefix}/genome_db/checkv_reps.tsv
    touch ${prefix}/genome_db/checkv_reps.dmnd
    mkdir -p ${prefix}/hmm_db/
    touch ${prefix}/hmm_db/test.hmm
    touch ${prefix}/hmm_db/checkv_hmms.tsv
    touch ${prefix}/hmm_db/genome_lengths.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        checkv: \$(checkv -h 2>&1  | sed -n 's/^.*CheckV v//; s/: assessing.*//; 1p')
    END_VERSIONS
    """
}
