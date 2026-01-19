process DOUBLET_DETECTION {
    tag "$meta.id"
    label 'process_medium'

    container "community.wave.seqera.io/library/scanpy_scrublet:latest"

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("*_doublets.h5ad"), emit: h5ad
    path "*.{png,pdf}"                       , emit: plots
    path "doublet_scores.csv"                , emit: scores
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    doublet_detection.py \\
        --input ${h5ad} \\
        --output ${prefix}_doublets.h5ad \\
        --scores doublet_scores.csv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
        scrublet: \$(python -c "import scrublet; print(scrublet.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_doublets.h5ad
    touch doublet_histogram.png
    touch doublet_umap.png
    touch doublet_scores.csv
    touch versions.yml
    """
}
