process LOGAN_DOWNLOAD {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/awscli:1.32.101--263008d47506de8c' :
        'community.wave.seqera.io/library/awscli:1.32.101--bca0ff06f80a537e' }"

    input:
    tuple val(meta), path(accessions)

    output:
    tuple val(meta), path("${prefix}.fa.zst")   , emit: fasta
    path "versions.yml"                         , emit: versions

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir tmp
    cat ${accessions} | xargs -I{} -n 1 -P ${task.cpus} aws s3 cp s3://logan-pub/c/{}/{}.contigs.fa.zst tmp --no-sign-request
    cat tmp/*.fa.zst > ${prefix}.fa.zst
    rm -rf tmp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awscli: \$( aws --version | sed 's/aws-cli\\///; s/ Python.*//' )
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fa.zst

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awscli: \$( aws --version | sed 's/aws-cli\\///; s/ Python.*//' )
    END_VERSIONS
    """
}
