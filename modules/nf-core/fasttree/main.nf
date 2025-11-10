process FASTTREE {
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fasttree:2.2.0--h7b50bb2_0' :
        'biocontainers/fasttree:2.2.0--h7b50bb2_0' }"

    input:
    tuple val(meta), path(alignment)

    output:
    tuple val(meta), path("*.tree"), emit: phylogeny
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args         = task.ext.args   ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def alignment    = alignment       ? "-nt ${alignment})" : ''
    """
    fasttree \\
        $args \\
        -log ${prefix}.tree.log \\
        ${alignment} \\
        > ${prefix}.tree

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fasttree: \$(fasttree -help 2>&1 | head -1  | sed 's/^FastTree \\([0-9\\.]*\\) .*\$/\\1/')
    END_VERSIONS
    """

    stub:
    def args         = task.ext.args   ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def alignment    = alignment       ? "-nt ${alignment})" : ''
    """
    touch ${prefix}.tree

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fasttree: \$(fasttree -help 2>&1 | head -1  | sed 's/^FastTree \\([0-9\\.]*\\) .*\$/\\1/')
    END_VERSIONS
    """
}
