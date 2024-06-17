process FILTERSEQUENCES {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-80c23cbcd32e2891421c54d1899665046feb07ef:77a31e289d22068839533bf21f8c4248ad274b60-0' :
        'biocontainers/mulled-v2-80c23cbcd32e2891421c54d1899665046feb07ef:77a31e289d22068839533bf21f8c4248ad274b60-0' }"

    input:
    tuple val(meta)     , path(fasta)
    tuple val(meta2)    , path(genomad_summary)
    tuple val(meta3)    , path(quality_summary)
    tuple val(meta4)    , path(tantan)
    tuple val(meta5)    , path(nuc_stats)
    path(contigs_to_keep)

    output:
    tuple val(meta), path("${prefix}_filtered.fasta.gz")    , emit: fasta
    tuple val(meta), path("${prefix}_filtering_data.tsv")   , emit: tsv
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def contigs_to_keep = contigs_to_keep ? "--contigs_to_keep ${contigs_to_keep}" : "--contigs_to_keep NO_CONTIGS_TO_KEEP_FILE"
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    if [ \$(zcat ${fasta} | grep ">" | wc -l) -gt 0 ]; then
        filtersequences.py \\
            --input ${fasta} \\
            --genomad_summary ${genomad_summary} \\
            --quality_summary ${quality_summary} \\
            --tantan ${tantan} \\
            --nucleotide_stats ${nuc_stats} \\
            --output_fasta ${prefix}_filtered.fasta \\
            --output_tsv ${prefix}_filtering_data.tsv \\
            ${contigs_to_keep} \\
            ${args}

        gzip ${prefix}_filtered.fasta
    else
        echo "" | gzip > ${prefix}_filtered.fasta.gz
        touch ${prefix}_filtering_data.tsv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed 's/Python //' )
        biopython: \$(echo \$(biopython_version.py 2>&1))
        pandas: \$(echo \$(pandas_version.py 2>&1))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo "" | gzip > ${prefix}_filtered.fasta.gz
    touch ${prefix}_filtering_data.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$( python --version | sed 's/Python //' )
        biopython: \$(echo \$(biopython_version.py 2>&1))
        pandas: \$(echo \$(pandas_version.py 2>&1))
    END_VERSIONS
    """
}
