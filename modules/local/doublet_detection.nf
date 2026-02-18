process DOUBLET_DETECTION {
    tag "$meta.id"
    label 'process_medium'

    conda "\"scanpy>=1.9\" \"scrublet>=0.2\""
    container "oras://community.wave.seqera.io/library/scrublet_scanpy:31ec5c9d8d3579e3"

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("*_doublets.h5ad"), emit: h5ad
    path "*.png"                             , emit: plots
    path "doublet_scores.csv"                , emit: scores
    path "doublet_attrition.json"            , emit: attrition, optional: true
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    // SLC: Build doublet detection parameters
    def remove_doublets = params.doublet_removal ? '--remove-doublets' : ''
    def save_attrition = params.save_attrition_log ? '--save-attrition-log' : ''
    def doublet_rate = params.expected_doublet_rate ?: 0.05
    
    """
    doublet_detection.py \\
        --input ${h5ad} \\
        --output ${prefix}_doublets.h5ad \\
        --scores doublet_scores.csv \\
        --expected-doublet-rate ${doublet_rate} \\
        ${remove_doublets} \\
        ${save_attrition} \\
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
    touch doublet_attrition.json
    touch versions.yml
    """
}
