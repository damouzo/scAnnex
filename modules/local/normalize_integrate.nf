process NORMALIZE_INTEGRATE {
    tag "$meta.id"
    label 'process_high'

    conda "\"scanpy>=1.9\" \"harmonypy>=0.0.9\""
    container "oras://community.wave.seqera.io/library/harmonypy_scanpy:a8efcfdf23c8acc8"

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("*_normalized.h5ad"), emit: h5ad
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    normalize_integrate.py \\
        --input ${h5ad} \\
        --output ${prefix}_normalized.h5ad \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
        harmonypy: \$(python -c "import harmonypy; print(harmonypy.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_normalized.h5ad
    touch versions.yml
    """
}
