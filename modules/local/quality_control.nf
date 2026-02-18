process QUALITY_CONTROL {
    tag "$meta.id"
    label 'process_medium'

    conda "\"scanpy>=1.9\""
    container "oras://community.wave.seqera.io/library/pip_scanpy:46ad0720691ef95a"

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("*_qc.h5ad")   , emit: h5ad
    path "qc_results/"                   , emit: qc_dir
    path "qc_results/*.png"              , emit: plots
    path "qc_results/qc_report.json"     , emit: report
    path "qc_results/cell_attrition_log.csv", emit: attrition_log, optional: true
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    // SLC: Build filtering arguments
    def filtering_method = ''
    if (params.use_quantile_filtering) {
        filtering_method = """--use-quantile-filtering \\
        --feature-quantile-low ${params.feature_quantile_low} \\
        --feature-quantile-high ${params.feature_quantile_high} \\
        --count-quantile-low ${params.count_quantile_low} \\
        --count-quantile-high ${params.count_quantile_high}"""
    } else if (params.use_mad_thresholds) {
        filtering_method = """--use-mad-thresholds \\
        --mad-threshold ${params.mad_multiplier}"""
    }
    
    def attrition_log = params.save_attrition_log ? '--save-attrition-log' : ''
    def max_mito = params.max_mito_percent ? "--max-mito ${params.max_mito_percent}" : ''
    
    """
    quality_control.py \\
        --input ${h5ad} \\
        --output ${prefix}_qc.h5ad \\
        --qc-dir qc_results \\
        --min-genes ${params.min_genes} \\
        --min-cells ${params.min_cells} \\
        ${max_mito} \\
        ${filtering_method} \\
        ${attrition_log} \\
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
    mkdir -p qc_results
    touch ${prefix}_qc.h5ad
    touch qc_results/qc_before_violin.png
    touch qc_results/qc_after_violin.png
    touch qc_results/qc_report.json
    touch qc_results/cell_attrition_log.csv
    touch versions.yml
    """
}
