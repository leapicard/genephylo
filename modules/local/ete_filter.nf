process ETE_FILTER {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::ete3==3.1.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ete3:3.1.2' :
        'quay.io/biocontainers/ete3:3.1.2' }"

    input:
    tuple val(meta), path(fasta)
    tuple val(meta2), path(taxidmap)

    output:
    tuple val(meta), path("*_filtered.fasta"), emit: fasta
    tuple val(meta2), path("*_speciescodes.txt"), emit: txt
    path "versions.yml", emit: versions

    when:
    taxdb_ready

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

    rename_seqs.py --taxidmap "$taxidmap" --input "$fasta" --prefix "${prefix}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ete3: \$(python -c "import ete3; print(ete3.__version__)")
    END_VERSIONS
    """
}