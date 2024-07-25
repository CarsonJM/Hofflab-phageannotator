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
    tuple val(meta), path("${prefix}.contigs.fa.gz")    , emit: contigs
    tuple val(meta), path("${prefix}.graph.fastg.gz")   , emit: graph
    tuple val(meta), path("${prefix}.megahit.log")      , emit: log
    path "versions.yml"                                 , emit: versions
    env min_kmer                                        , emit: min_kmer
    env max_kmer                                        , emit: max_kmer
    // tuple val(meta), path("megahit_out/intermediate_contigs/k*.contigs.fa.gz")      , emit: k_contigs
    // tuple val(meta), path("megahit_out/intermediate_contigs/k*.addi.fa.gz")         , emit: addi_contigs
    // tuple val(meta), path("megahit_out/intermediate_contigs/k*.local.fa.gz")        , emit: local_contigs
    // tuple val(meta), path("megahit_out/intermediate_contigs/k*.final.contigs.fa.gz"), emit: kfinal_contigs
    

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
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
        megahit_toolkit contig2fastg \$kmer_size megahit_out/${prefix}.contigs.fa > ${prefix}.graph.fastg

        mv megahit_out/*.log ${prefix}.megahit.log
        gzip -c megahit_out/*.contigs.fa > ${prefix}.contigs.fa.gz
        gzip ${prefix}.graph.fastg

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
        megahit_toolkit contig2fastg \$kmer_size megahit_out/${prefix}.contigs.fa > ${prefix}.graph.fastg

        mv megahit_out/*.log ${prefix}.megahit.log
        gzip -c megahit_out/*.contigs.fa > ${prefix}.contigs.fa.gz
        gzip ${prefix}.graph.fastg

        # identify min/max kmer size
        kmer_string=\$(grep "k list: " ${prefix}.megahit.log | sed 's/.*k list: //; s/ .*//')
        kmer_array=(\${kmer_string//,/ })
        min_kmer=\$(IFS=\$'\\n'; echo "\${kmer_array[*]}" | sort -nr | tail -n1)
        max_kmer=\$(IFS=\$'\\n'; echo "\${kmer_array[*]}" | sort -nr | head -n1)

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
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.megahit.log
    echo "" | gzip -c > ${prefix}.contigs.fa.gz
    echo "" | gzip -c > ${prefix}.graph.fastg.gz
    min_kmer=21
    max_kmer=141

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        megahit: \$(echo \$(megahit -v 2>&1) | sed 's/MEGAHIT v//')
    END_VERSIONS
    """
}
