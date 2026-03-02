process MERGE_SAMPLES {
    tag "merge_${h5ad_list.size()}_samples"
    label 'process_medium'

    conda "\"scanpy>=1.9\" \"scipy\""
    container "oras://community.wave.seqera.io/library/scanpy_scipy:af35be00f10024f0"

    input:
    path(h5ad_list)

    output:
    path("merged_samples.h5ad"), emit: h5ad
    path "versions.yml"        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def h5ad_files = h5ad_list.collect { it.toString() }.join(' ')
    """
    merge_samples.py \\
        --inputs ${h5ad_files} \\
        --output merged_samples.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
        scipy: \$(python -c "import scipy; print(scipy.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch merged_samples.h5ad
    touch versions.yml
    """
}
