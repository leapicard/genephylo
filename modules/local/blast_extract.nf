process BLAST_EXTRACT {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::gawk=5.3.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.1' :
        'quay.io/biocontainers/gawk:5.3.1' }"

    input:
    tuple val(meta), path(blastres)

    output:
    tuple val(meta), path("*_accessions.txt")       , emit: accessions
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gawk '{print \$1}' "$blastres" > "${prefix}_accessions.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gawk: \$(gawk --version | head -n 1 | cut -d',' -f1 | cut -d' ' -f3)
    END_VERSIONS
    """
}