process PREPARE_DASHBOARD {
    label 'process_low'
    
    container 'community.wave.seqera.io/library/scanpy_pandas:latest'

    input:
    path(h5ad_files)

    output:
    path "dashboard/*", emit: dashboard
    path "versions.yml", emit: versions

    script:
    """
    mkdir -p dashboard

    prepare_dashboard.py \\
        --input ${h5ad_files} \\
        --output_dir dashboard

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p dashboard
    touch dashboard/metadata.csv
    touch dashboard/umap_coords.csv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
    END_VERSIONS
    """
}
