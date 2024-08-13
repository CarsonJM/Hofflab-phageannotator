process COBRAMETA {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/cobra-meta:1.2.3--10daf775a65e0312':
        'biocontainers/cobra-meta:1.2.3--pyhdfd78af_0' }"

    input:
    tuple val(meta) , path(fasta)
    tuple val(meta2), path(coverage)
    tuple val(meta3), path(query)
    tuple val(meta4), path(bam)
    val assembler
    val mink
    val maxk

    output:
    tuple val(meta), path("${prefix}_COBRA_extended.fasta.gz")      , emit: fasta
    tuple val(meta), path("${prefix}_COBRA_joining_summary.txt")    , emit: joining_summary
    tuple val(meta), path("${prefix}_log")                          , emit: log
    path "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    prefix      = task.ext.prefix ?: "${meta.id}"
    fasta_name  = fasta.getName().replace(".gz", "")
    """
    gunzip -c ${fasta} > ${fasta_name}
    tail ${virus_summary} -n +2 | awk '{print \$1}' > ${prefix}_viral_assemblies.txt
    tail ${coverage} -n +2 > ${prefix}_cobra_coverage.txt

    cobra-meta \\
        --fasta ${fasta_name} \\
        --coverage ${prefix}_cobra_coverage.txt \\
        --query ${prefix}_viral_assemblies.txt \\
        --mapping ${bam} \\
        --assembler ${assembler} \\
        --mink ${mink} \\
        --maxk ${maxk} \\
        --threads ${task.cpus} \\
        --output ${prefix} \\
        $args

    if [ -f ${prefix}/COBRA_extended.fasta ]; then
        mv ${prefix}/COBRA_extended.fasta ${prefix}_COBRA_extended.fasta
        gzip ${prefix}_COBRA_extended.fasta
        mv ${prefix}/COBRA_joining_summary.txt ${prefix}_COBRA_joining_summary.txt
        mv ${prefix}/log ${prefix}_log
    else 
        echo "" | gzip > ${prefix}_COBRA_extended.fasta.gz
        touch ${prefix}_COBRA_joining_summary.txt
        mv ${prefix}/log ${prefix}_log
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cobra: \$(echo \$(cobra-meta --version 2>&1) | sed 's/^.*cobra v//' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}_COBRA_extended.fasta.gz
    touch ${prefix}_COBRA_joining_summary.txt
    touch ${prefix}_log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cobra: \$(echo \$(cobra-meta --version 2>&1) | sed 's/^.*cobra v//' ))
    END_VERSIONS
    """
}
