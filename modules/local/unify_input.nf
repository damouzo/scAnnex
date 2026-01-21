process UNIFY_INPUT {
    tag "$meta.id"
    label 'process_low'

    conda "bioconda::scanpy=1.10.0 bioconda::anndata=0.10.3 conda-forge::numpy=1.24.0"
    container "quay.io/biocontainers/scanpy:1.10.0--pyhdfd78af_0"

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("*.h5ad"), emit: h5ad
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sample_id = meta.id
    def batch = meta.batch ?: 'batch_unknown'
    def condition = meta.condition ?: 'condition_unknown'
    def input_type = meta.file_type ?: 'h5ad'
    
    """
    unify_input.py \\
        --input ${input_file} \\
        --input-type ${input_type} \\
        --output ${prefix}_unified.h5ad \\
        --sample-id ${sample_id} \\
        --batch ${batch} \\
        --condition ${condition} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
        anndata: \$(python -c "import anndata; print(anndata.__version__)")
        python: \$(python --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_unified.h5ad
    touch versions.yml
    """
}
