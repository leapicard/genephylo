process ETE_TAXDB {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::ete3==3.1.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ete3:3.1.2' :
        'quay.io/biocontainers/ete3:3.1.2' }"

    output:
    path "taxdb_update.done"
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Set up writable environment for ete3
    export HOME=\$PWD
    export XDG_DATA_HOME=\$PWD/.local/share
    export XDG_CONFIG_HOME=\$PWD/.config
    export XDG_CACHE_HOME=\$PWD/.cache
    
    # Create necessary directories
    mkdir -p .local/share .config .cache .etetoolkit

    python -c "from ete3 import NCBITaxa; NCBITaxa().update_taxonomy_database()"
    echo "ETE3 taxonomy database updated on $(date -u)" > taxdb_update.done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ete3: \$(python -c "import ete3; print(ete3.__version__)")
    END_VERSIONS
    """
}