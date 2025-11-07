process BLAST_FILTER {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::biopython=1.76"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.76--1' :
        'quay.io/biocontainers/biopython:1.76--1' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(taxidmap)

    output:
    tuple val(meta), path("*_renamed.fasta"), emit: fasta
    tuple val(meta2), path("*_speciescodes.txt"), emit: txt
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    beforeScript "export XDG_CONFIG_HOME=${task.workDir}/.config"

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    rename_seqs.py --taxidmap "$taxidmap" --input "$fasta" --prefix "${prefix}" --threads ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | cut -d' ' -f2)
    END_VERSIONS
    """
}