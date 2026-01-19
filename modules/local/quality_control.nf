process QUALITY_CONTROL {
    tag "$meta.id"
    label 'process_medium'

    container "community.wave.seqera.io/library/scanpy_anndata:latest"

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("*_qc.h5ad")   , emit: h5ad
    path "*.{png,pdf}"                   , emit: plots
    path "qc_metrics.csv"                , emit: metrics
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    quality_control.py \\
        --input ${h5ad} \\
        --output ${prefix}_qc.h5ad \\
        --metrics qc_metrics.csv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_qc.h5ad
    touch qc_violin.png
    touch qc_scatter.png
    touch qc_metrics.csv
    touch versions.yml
    """
}
