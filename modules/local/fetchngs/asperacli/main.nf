process ASPERACLI {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/aspera-cli:4.14.0--hdfd78af_1' :
        'biocontainers/aspera-cli:4.14.0--hdfd78af_1' }"

    input:
    tuple val(meta), val(fastq)
    val user

    output:
    tuple val(meta), path("*fastq.gz"), emit: fastq
    tuple val(meta), path("*md5")     , emit: md5
    path "versions.yml"               , emit: versions

    script:
    def args = task.ext.args ?: ''
    def conda_prefix = ['singularity', 'apptainer'].contains(workflow.containerEngine) ? "export CONDA_PREFIX=/usr/local" : ""
    prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        $conda_prefix

        ascp \\
            $args \\
            -i \$CONDA_PREFIX/etc/aspera/aspera_bypass_dsa.pem \\
            ${user}@${fastq[0]} \\
            ${prefix}.fastq.gz

        echo "${meta.md5_1}  ${prefix}.fastq.gz" > ${prefix}.fastq.gz.md5
        md5sum -c ${prefix}.fastq.gz.md5

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            aspera_cli: \$(ascli --version)
        END_VERSIONS
        """
    } else {
        """
        $conda_prefix

        ascp \\
            $args \\
            -i \$CONDA_PREFIX/etc/aspera/aspera_bypass_dsa.pem \\
            ${user}@${fastq[0]} \\
            ${prefix}_1.fastq.gz

        echo "${meta.md5_1}  ${prefix}_1.fastq.gz" > ${prefix}_1.fastq.gz.md5
        md5sum -c ${prefix}_1.fastq.gz.md5

        ascp \\
            $args \\
            -i \$CONDA_PREFIX/etc/aspera/aspera_bypass_dsa.pem \\
            ${user}@${fastq[1]} \\
            ${prefix}_2.fastq.gz

        echo "${meta.md5_2}  ${prefix}_2.fastq.gz" > ${prefix}_2.fastq.gz.md5
        md5sum -c ${prefix}_2.fastq.gz.md5

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            aspera_cli: \$(ascli --version)
        END_VERSIONS
        """
    }


    stub:
    def args = task.ext.args ?: ''
    def conda_prefix = ['singularity', 'apptainer'].contains(workflow.containerEngine) ? "export CONDA_PREFIX=/usr/local" : ""
    prefix = task.ext.prefix ?: "${meta.id}"
    if (meta.single_end) {
        """
        echo "" | gzip > ${prefix}.fastq.gz
        touch ${prefix}.fastq.gz.md5

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            aspera_cli: \$(ascli --version)
        END_VERSIONS
        """
    } else {
        """
        echo "" | gzip > ${prefix}_1.fastq.gz
        echo "" | gzip > ${prefix}_2.fastq.gz
        touch ${prefix}.fastq.gz.md5

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            aspera_cli: \$(ascli --version)
        END_VERSIONS
        """
    }
}
