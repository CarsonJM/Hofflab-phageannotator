process PHABLES_RUN {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://carsonjm/phables:1.4.0':
        'docker://carsonjm/phables:1.4.0' }"

    input:
    tuple val(meta) , path(graph)
    tuple val(meta2), path(fastq)
    path phables_db

    output:
    tuple val(meta), path("${prefix}_phables.fasta.gz") , emit: fasta
    tuple val(meta), path("${prefix}_phables.log")      , emit: log
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    bash /usr/local/bin/_entrypoint.sh
    mkdir reads
    cp ${fastq[0]} reads/
    cp ${fastq[1]} reads/
    mkdir phables.out

    phables \\
        run \\
        --configfile ${projectDir}/assets/phables/config.yml \\
        --input ${graph} \\
        --reads reads \\
        --threads ${task.cpus} \\
        --prefix ${prefix} \\
        ${args}

    mv phables_out/phables/resolved_paths.fasta ./${prefix}_phables.fasta
    gzip ./${prefix}_phables.fasta
    mv phables_out/phables.log ./${prefix}_phables.log
    rm -rf phables_out/

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phables: \$(echo \$(phables -v 2>&1) | sed -n 's/phables, version //' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_phables.log
    echo "" | gzip > ${prefix}_phables.fasta.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phables: \$(echo \$(phables -v 2>&1) | sed -n 's/phables, version //' ))
    END_VERSIONS
    """
}
