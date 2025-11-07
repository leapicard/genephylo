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

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Set up writable environment for ete3
    export HOME=\$PWD
    export XDG_DATA_HOME=\$PWD/.local/share
    export XDG_CONFIG_HOME=\$PWD/.config
    export XDG_CACHE_HOME=\$PWD/.cache
    
    # Create necessary directories
    mkdir -p .local/share .config .cache .etetoolkit

    rename_seqs.py --taxidmap "$taxidmap" --input "$fasta" --prefix "${prefix}" --workdir \$PWD

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ete3: \$(ete3 version | cut -d' ' -f1)
    END_VERSIONS
    """
}