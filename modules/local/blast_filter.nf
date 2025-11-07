process BLAST_FILTER {
    tag "$meta.id"
    label 'process_medium'

    conda "conda-forge::ete3==3.1.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ete3:3.1.2' :
        'quay.io/biocontainers/ete3:3.1.2' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(taxidmap)

    output:
    tuple val(meta), path("*_renamed.fasta"), emit: fasta
    tuple val(meta2), path("*_speciescodes.txt"), emit: txt
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    env:
    XDG_CACHE_HOME = "\$PWD/.cache"
    HOME = "\$PWD"

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    rename_seqs.py --taxidmap "$taxidmap" --input "$fasta" --prefix "${prefix}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ete3: \$(ete3 version | cut -d' ' -f1)
    END_VERSIONS
    """
}