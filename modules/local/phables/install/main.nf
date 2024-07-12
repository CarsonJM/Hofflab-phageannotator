process PHABLES_INSTALL {
    tag "phables_install"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/phables:1.4.0--pyhdfd78af_0':
        'biocontainers/phables:1.4.0--pyhdfd78af_0' }"

    output:
    path "phables_db"           , emit: phables_db
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    mkdir phables_db

    wget -P phables_db https://raw.githubusercontent.com/metagentools/MetaCoAG/develop/src/metacoag/metacoag_utils/auxiliary/marker.hmm
    wget -P phables_db https://phrogs.lmge.uca.fr/downloads_from_website/phrog_annot_v4.tsv

    wget -P phables_db https://phrogs.lmge.uca.fr/downloads_from_website/phrogs_mmseqs_db.tar.gz
    tar -xvf phables_db/phrogs_mmseqs_db.tar.gz -C phables_db


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phables: \$(echo \$(phables -v 2>&1) | sed -n 's/phables, version //' ))
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    """
    mkdir phables_db
    touch phables_db/marker.hmm
    touch phables_db/phrog_annot_v4.tsv
    touch phables_db/phrogs_mmseqs_db.tar.gz


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        phables: \$(echo \$(phables -v 2>&1) | sed -n 's/phables, version //' ))
    END_VERSIONS
    """
}