process MEGAHIT {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/megahit_pigz:657d77006ae5f222' :
        'community.wave.seqera.io/library/megahit_pigz:87a590163e594224' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}.megahit.fa.gz")    , emit: contigs
    tuple val(meta), path("${prefix}.megahit.fastg.gz") , emit: graph
    tuple val(meta), path("${prefix}.megahit.log")      , emit: log
    tuple val(meta), env(min_kmer)                      , emit: min_kmer
    tuple val(meta), env(max_kmer)                      , emit: max_kmer
    path "versions.yml"                                 , emit: versions
    // tuple val(meta), path("intermediate_contigs/k*.contigs.fa.gz")      , emit: k_contigs
    // tuple val(meta), path("intermediate_contigs/k*.addi.fa.gz")         , emit: addi_contigs
    // tuple val(meta), path("intermediate_contigs/k*.local.fa.gz")        , emit: local_contigs
    // tuple val(meta), path("intermediate_contigs/k*.final.contigs.fa.gz"), emit: kfinal_contigs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def reads_command = meta.single_end ? "-r ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    megahit \\
        ${reads_command} \\
        ${args} \\
        -t ${task.cpus} \\
        --out-prefix ${prefix}

    # create assembly graph file
    kmer_size=\$(grep "^>" megahit_out/${prefix}.contigs.fa | sed 's/>k//; s/_.*//')
    megahit_toolkit contig2fastg \$kmer_size megahit_out/${prefix}.contigs.fa > ${prefix}.megahit.fastg

    # rename output files
    mv megahit_out/*.log ${prefix}.megahit.log
    gzip -c megahit_out/*.contigs.fa > ${prefix}.megahit.fa.gz
    gzip ${prefix}.megahit.fastg

    # identify min/max kmer size
    kmer_string=\$(grep "k list: " ${prefix}.megahit.log | sed 's/.*k list: //; s/ .*//')
    kmer_array=(\${kmer_string//,/ })
    min_kmer=\$(IFS=\$'\\n'; echo "\${kmer_array[*]}" | sort -nr | tail -n 1)
    max_kmer=\$(IFS=\$'\\n'; echo "\${kmer_array[*]}" | sort -nr | head -n 1)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        megahit: \$(echo \$(megahit -v 2>&1) | sed 's/MEGAHIT v//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def reads_command = meta.single_end ? "-r ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p intermediate_contigs
    echo "" | gzip > ${prefix}.megahit.fa.gz
    echo "" | gzip > ${prefix}.megahit.fastg.gz
    touch ${prefix}.megahit.log
    min_kmer=21
    max_kmer=141

    #echo "" | gzip > intermediate_contigs/k21.contigs.fa.gz
    #echo "" | gzip > intermediate_contigs/k21.addi.fa.gz
    #echo "" | gzip > intermediate_contigs/k21.local.fa.gz
    #echo "" | gzip > intermediate_contigs/k21.final.contigs.fa.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        megahit: \$(echo \$(megahit -v 2>&1) | sed 's/MEGAHIT v//')
    END_VERSIONS
    """
}
