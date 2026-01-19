process DIMENSIONALITY_REDUCTION {
    tag "$meta.id"
    label 'process_high'

    container "community.wave.seqera.io/library/scanpy_leidenalg:latest"

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("*_processed.h5ad"), emit: h5ad
    path "*.{png,pdf}"                        , emit: plots
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    dimensionality_reduction.py \\
        --input ${h5ad} \\
        --output ${prefix}_processed.h5ad \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
        leidenalg: \$(python -c "import leidenalg; print(leidenalg.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_processed.h5ad
    touch pca_variance.png
    touch umap_clusters.png
    touch versions.yml
    """
}
