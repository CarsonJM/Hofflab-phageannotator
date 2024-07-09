process IPHOP_PREDICT {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '/mmfs1/gscratch/pedslabs_hoffman/carsonjm/iphop_1.3.3--pyhdfd78af_0':
        '/mmfs1/gscratch/pedslabs_hoffman/carsonjm/iphop_1.3.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path iphop_db

    output:
    tuple val(meta), path("Host_prediction_to_genus_m*.csv")    , emit: iphop_genus
    tuple val(meta), path("Host_prediction_to_genome_m*.csv")   , emit: iphop_genome
    tuple val(meta), path("Detailed_output_by_tool.csv")        , emit: iphop_detailed_output
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def min_score = args.contains('--min_score') ? args.split('--min_score ')[1] : '90'
    prefix = task.ext.prefix ?: "${meta.id}"
    fasta_name  = fasta.getName().replace(".gz", "")
    """
    if [ \$(zcat ${fasta} | grep ">" | wc -l) -gt 0 ]; then
        gunzip -c ${fasta} > ${fasta_name}

        iphop \\
            predict \\
            --fa_file ${fasta_name} \\
            --out_dir iphop_results \\
            --db_dir ${iphop_db} \\
            --num_threads ${task.cpus} \\
            $args

        mv iphop_results/Host_prediction_to_genus_m*.csv .
        mv iphop_results/Host_prediction_to_genome_m*.csv .
        mv iphop_results/Detailed_output_by_tool.csv .
    else
        mkdir -p iphop_results
        touch Host_prediction_to_genus_m${min_score}.csv
        touch Host_prediction_to_genome_m${min_score}.csv
        touch Detailed_output_by_tool.csv
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iphop: \$(echo \$(iphop --version 2>&1) | head -n 1 | sed 's/^.*iPHoP v//; s/: integrating.*\$//' )
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def min_score = args.contains('--min_score') ? args.split('--min_score ')[1] : '90'
    """
    mkdir -p iphop_results
    touch Host_prediction_to_genus_m${min_score}.csv
    touch Host_prediction_to_genome_m${min_score}.csv
    touch Detailed_output_by_tool.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        iphop: \$(echo \$(iphop --version 2>&1) | head -n 1 | sed 's/^.*iPHoP v//; s/: integrating.*\$//' )
    END_VERSIONS
    """
}
