process GENOMAD_ENDTOEND {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/genomad:1.8.0--pyhdfd78af_0':
        'biocontainers/genomad:1.8.0--pyhdfd78af_0' }"

    input:
    tuple val(meta) , path(fasta)
    path(genomad_db)

    output:
    tuple val(meta), path("*_virus_summary.tsv")                , emit: virus_summary
    tuple val(meta), path("*_virus.fna.gz")                     , emit: virus_fasta
    // tuple val(meta), path("*_aggregated_classification.tsv")    , emit: aggregated_classification
    // tuple val(meta), path("*_taxonomy.tsv")                     , emit: taxonomy
    // tuple val(meta), path("*_genes.tsv")                        , emit: genes
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    filename = fasta.toString().tokenize('.')[0]
    """
    if [ \$(zcat ${fasta} | grep ">" | wc -l) -gt 0 ]; then
        genomad \\
            end-to-end \\
            ${fasta} \\
            ./ \\
            ${genomad_db} \\
            --threads ${task.cpus} \\
            ${args}

        mv *_aggregated_classification/*_aggregated_classification.tsv .
        mv *_annotate/*_taxonomy.tsv .
        mv *_annotate/*_genes.tsv .
        mv *_summary/*_virus_summary.tsv .
        mv *_summary/*_virus.fna .
        gzip *_virus.fna

        rm -rf ./${filename}_*/*
    else
        touch ${filename}_aggregated_classification.tsv
        touch ${filename}_taxonomy.tsv
        touch ${filename}_genes.tsv
        touch ${filename}_virus_summary.tsv
        echo "" | gzip > ${filename}_virus.fna.gz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genomad: \$(echo \$(genomad --version 2>&1) | sed 's/^.*geNomad, version //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def filename = "${fasta}"[0..<"${fasta}".lastIndexOf('.')]
    """
    touch ${filename}_aggregated_classification.tsv
    touch ${filename}_taxonomy.tsv
    touch ${filename}_genes.tsv
    touch ${filename}_virus_summary.tsv
    echo "" | gzip > ${filename}_virus.fna.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genomad: \$(echo \$(genomad --version 2>&1) | sed 's/^.*geNomad, version //; s/ .*\$//')
    END_VERSIONS
    """
}
