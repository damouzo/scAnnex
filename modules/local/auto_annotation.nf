process AUTO_ANNOTATION {
    tag "$meta.id"
    label 'process_medium'

    container "community.wave.seqera.io/library/scanpy_pandas:latest"

    input:
    tuple val(meta), path(h5ad)
    path marker_list

    output:
    tuple val(meta), path("*_annotated.h5ad"), emit: h5ad
    path "*.{png,pdf}"                        , emit: plots
    path "annotation_summary.csv"             , emit: summary
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    auto_annotation.py \\
        --input ${h5ad} \\
        --marker-list ${marker_list} \\
        --output ${prefix}_annotated.h5ad \\
        --summary annotation_summary.csv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_annotated.h5ad
    touch dotplot_markers.png
    touch umap_celltypes.png
    touch annotation_summary.csv
    touch versions.yml
    """
}
