process MEGAHIT {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-0f92c152b180c7cd39d9b0e6822f8c89ccb59c99:8ec213d21e5d03f9db54898a2baeaf8ec729b447-0' :
        'biocontainers/mulled-v2-0f92c152b180c7cd39d9b0e6822f8c89ccb59c99:8ec213d21e5d03f9db54898a2baeaf8ec729b447-0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("megahit_out/*.contigs.fa.gz")                            , emit: contigs
    tuple val(meta), path("megahit_out/*.graph.fastg.gz")                           , emit: graph
    tuple val(meta), path("megahit_out/intermediate_contigs/k*.contigs.fa.gz")      , emit: k_contigs
    tuple val(meta), path("megahit_out/intermediate_contigs/k*.addi.fa.gz")         , emit: addi_contigs
    tuple val(meta), path("megahit_out/intermediate_contigs/k*.local.fa.gz")        , emit: local_contigs
    tuple val(meta), path("megahit_out/intermediate_contigs/k*.final.contigs.fa.gz"), emit: kfinal_contigs
    path "versions.yml"                                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def readList = reads instanceof List ? reads.collect{ it.toString() } : [reads.toString()]
    if ( meta.single_end ) {
        """
        megahit \\
            -r ${readList.join(',')} \\
            -t $task.cpus \\
            $args \\
            --out-prefix $prefix

        # create assembly graph file
        kmer_size=\$(grep "^>" megahit_out/${prefix}.contigs.fa | sed 's/>k//; s/_.*//')
        megahit_toolkit contig2fastg \$kmer_size megahit_out/${prefix}.contigs.fa > megahit_out/${prefix}.graph.fastg

        pigz \\
            --no-name \\
            -p $task.cpus \\
            $args2 \\
            megahit_out/*.fa \\
            megahit_out/*.fastg \\
            megahit_out/intermediate_contigs/*.fa

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            megahit: \$(echo \$(megahit -v 2>&1) | sed 's/MEGAHIT v//')
        END_VERSIONS
        """
    } else {
        def read1 = []
        def read2 = []
        readList.eachWithIndex{ v, ix -> ( ix & 1 ? read2 : read1 ) << v }
        """
        megahit \\
            -1 ${read1.join(',')} \\
            -2 ${read1.join(',')} \\
            -t $task.cpus \\
            $args \\
            --out-prefix $prefix
        
        # create assembly graph file
        kmer_size=\$(grep "^>" megahit_out/${prefix}.contigs.fa | sed 's/>k//; s/_.*//')
        megahit_toolkit contig2fastg \$kmer_size megahit_out/${prefix}.contigs.fa > megahit_out/${prefix}.graph.fastg

        pigz \\
            --no-name \\
            -p $task.cpus \\
            $args2 \\
            megahit_out/*.fa \\
            megahit_out/*.fastg \\
            megahit_out/intermediate_contigs/*.fa

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            megahit: \$(echo \$(megahit -v 2>&1) | sed 's/MEGAHIT v//')
        END_VERSIONS
        """
    }

    stub:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def reads_1 = reads[0].join(',')
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p megahit_out/intermediate_contigs
    echo "" | gzip > megahit_out/${prefix}.contigs.fa.gz
    echo "" | gzip > megahit_out/${prefix}.graph.fastg.gz
    echo "" | gzip > megahit_out/intermediate_contigs/ktest.contigs.fa.gz
    echo "" | gzip > megahit_out/intermediate_contigs/ktest.addi.fa.gz
    echo "" | gzip > megahit_out/intermediate_contigs/ktest.local.fa.gz
    echo "" | gzip > megahit_out/intermediate_contigs/ktest.final.contigs.fa.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        megahit: \$(echo \$(megahit -v 2>&1) | sed 's/MEGAHIT v//')
    END_VERSIONS
    """
}
