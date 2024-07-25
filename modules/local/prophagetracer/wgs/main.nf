process PROPHAGETRACER_WGS {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/blast_bwa_sambamba_samtools_gawk:6d2fecd74dcf1c21' :
        'community.wave.seqera.io/library/blast_bwa_sambamba_samtools:4b89bb5b2ebc9725' }"

    input:
    tuple val(meta) , path(fasta)
    tuple val(meta2), path(bam)

    output:
    tuple val(meta), path("${prefix}.prophage.out") , emit: prophage
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args   ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getExtension() == "gz" ? true : false
    def fasta_name = is_compressed ? fasta.getBaseName() : fasta
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    samtools view ${bam} -o ${prefix}.sam

    prophage_tracer_WGS.sh \\
        -m ${prefix}.sam \\
        -r ${fasta_name} \\
        -p ${prefix} \\
        -t ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prophagetracer: \$(echo \$( prophage_tracer_WGS.sh -h | sed -n "s/Prophage Tracer V//; s/\s.*//;  2p" ))
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.prophage.out

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prophagetracer: \$(echo \$( prophage_tracer_WGS.sh -h | sed -n "s/Prophage Tracer V//; s/\s.*//;  2p" ))
    END_VERSIONS
    """
}
