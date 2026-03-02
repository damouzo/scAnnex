process DIFFERENTIAL_EXPRESSION {
    tag "merged_samples"
    label 'process_medium'

    conda "\"scanpy>=1.9\" \"scipy\""
    container "oras://community.wave.seqera.io/library/scanpy_scipy:af35be00f10024f0"

    input:
    path h5ad
    path contrasts_file

    output:
    path "dge_results/"                      , emit: results_dir
    path "dge_results/*.csv"                 , emit: tables
    path "dge_results/plots/*.{png,pdf}"     , emit: plots, optional: true
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    // Build DGE arguments
    def contrasts_arg = contrasts_file.name != 'NO_FILE' ? "--contrasts ${contrasts_file}" : ''
    def groupby = params.dge_groupby ? "--groupby ${params.dge_groupby}" : ''
    def reference = params.dge_reference ? "--reference ${params.dge_reference}" : ''
    def method = params.dge_method ? "--method ${params.dge_method}" : '--method wilcoxon'
    def min_pct = params.dge_min_pct ? "--min-pct ${params.dge_min_pct}" : '--min-pct 0.1'
    def logfc_threshold = params.dge_logfc_threshold ? "--logfc-threshold ${params.dge_logfc_threshold}" : '--logfc-threshold 0.25'
    def pval_cutoff = params.dge_pval_cutoff ? "--pval-cutoff ${params.dge_pval_cutoff}" : '--pval-cutoff 0.05'
    def top_n = params.dge_top_n_genes ? "--top-n ${params.dge_top_n_genes}" : '--top-n 50'
    def save_plots = params.dge_save_plots ? '--save-plots' : ''
    
    """
    # Run DGE (output H5AD not required)
    differential_expression.py \\
        --input ${h5ad} \\
        --output merged_dge_results.h5ad \\
        --output-dir dge_results \\
        ${contrasts_arg} \\
        ${groupby} \\
        ${reference} \\
        ${method} \\
        ${min_pct} \\
        ${logfc_threshold} \\
        ${pval_cutoff} \\
        ${top_n} \\
        ${save_plots} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
        scipy: \$(python -c "import scipy; print(scipy.__version__)")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p dge_results/plots
    touch dge_results/dge_results_combined.csv
    touch dge_results/plots/volcano_plot.png
    touch versions.yml
    """
}
