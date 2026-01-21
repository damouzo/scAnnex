process STANDARD_PROCESSING {
    tag "$meta.id"
    label 'process_high'

    conda "bioconda::scanpy=1.10.0 bioconda::anndata=0.10.3 conda-forge::numpy=1.24.0"
    container "quay.io/biocontainers/scanpy:1.10.0--pyhdfd78af_0"

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("*_processed.h5ad"), emit: h5ad
    path "standard_processing_results/"        , emit: results_dir
    path "standard_processing_results/*.png"   , emit: plots
    path "standard_processing_results/umap_coordinates.csv", emit: umap_coords
    path "standard_processing_results/cell_metadata.csv"   , emit: metadata
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    // SLC: Build clustering parameters
    def resolutions = params.clustering_resolutions ?: '0.1,0.3,0.5,0.7,0.9'
    def default_res = params.default_clustering_resolution ?: 0.5
    def method = params.clustering_method ?: 'leiden'
    
    """
    standard_processing.py \\
        --input ${h5ad} \\
        --output ${prefix}_processed.h5ad \\
        --output-dir standard_processing_results \\
        --target-sum ${params.target_sum} \\
        --n-top-genes ${params.n_top_genes} \\
        --n-pcs ${params.n_pcs} \\
        --n-neighbors ${params.n_neighbors} \\
        --umap-min-dist ${params.umap_min_dist} \\
        --clustering-method ${method} \\
        --clustering-resolutions ${resolutions} \\
        --default-resolution ${default_res} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
        anndata: \$(python -c "import anndata; print(anndata.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p standard_processing_results
    touch ${prefix}_processed.h5ad
    touch standard_processing_results/pca_variance.png
    touch standard_processing_results/clustering_multi_resolution.png
    touch standard_processing_results/umap_coordinates.csv
    touch standard_processing_results/cell_metadata.csv
    touch versions.yml
    """
}
