process EXTRACT_ACCESSIONS {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::gawk=5.3.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0' :
        'quay.io/biocontainers/gawk:5.3.1' }"

    input:
    tuple val(meta), path(blastres)

    output:
    tuple val(meta), null, path("${prefix}_accessions.txt"), emit: accessions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    gawk '{print \$1}' "$blastres" > "${prefix}_accessions.txt"
    """
}