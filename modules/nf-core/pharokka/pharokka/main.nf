process PHAROKKA_PHAROKKA {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pharokka:1.5.1--pyhdfd78af_0':
        'biocontainers/pharokka:1.5.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path pharokka_db

    output:
    tuple val(meta), path("${prefix}_cds_final_merged_output.tsv")  , emit: cds_final_merged_output
    tuple val(meta), path("${prefix}_top_hits_mash_inphared.tsv")   , emit: mash                    , optional: true
    path "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    fasta_name  = fasta.getName().replace(".gz", "")
    """
    gunzip -c ${fasta} > ${fasta_name}

    pharokka.py \\
        --infile ${fasta_name} \\
        --outdir ${prefix}_pharokka \\
        --database ${pharokka_db} \\
        --threads ${task.cpus} \\
        --prefix ${prefix} \\
        $args
    mv ${prefix}_pharokka/${prefix}_top_hits_mash_inphared.tsv ${prefix}_top_hits_mash_inphared.tsv
    mv ${prefix}_pharokka/${prefix}_cds_final_merged_output.tsv ${prefix}_cds_final_merged_output.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pharokka: \$(pharokka.py --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"

    """
    mkdir -p ${prefix}_pharokka
    touch ${prefix}_pharokka/${prefix}.gbk
    touch ${prefix}_pharokka/${prefix}.log
    touch ${prefix}_pharokka/${prefix}_cds_functions.tsv
    touch ${prefix}_pharokka/${prefix}_top_hits_card.tsv
    touch ${prefix}_pharokka/top_hits_vfdb.tsv
    touch ${prefix}_pharokka/${prefix}_top_hits_inphared
    touch ${prefix}_pharokka/${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pharokka: \$(pharokka.py --version)
    END_VERSIONS
    """
}
