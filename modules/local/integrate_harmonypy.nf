process INTEGRATE_HARMONYPY {
    tag "integrate_${h5ad_list.size()}_samples"
    label 'process_high'

    conda "\"scanpy>=1.9\" \"harmonypy>=0.0.9\" \"scikit-learn\" \"scipy\""
    container "oras://community.wave.seqera.io/library/harmonypy_scanpy_scikit-learn_scipy:b98faf93ef84524e"

    input:
    path(h5ad_list)
    path split_script

    output:
    path "integrated.h5ad"               , emit: h5ad
    path "integrated_samples/*.h5ad"     , emit: sample_h5ad, optional: true
    path "integration_results/*"         , emit: integration_results
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def h5ad_files = h5ad_list.collect { h5ad_path -> h5ad_path.toString() }.join(' ')

    """
    merge_samples.py \\
        --inputs ${h5ad_files} \\
        --output merged_for_integration.h5ad

    normalize_integrate.py \\
        --input merged_for_integration.h5ad \\
        --output integrated.h5ad \\
        --run-integration \\
        --batch-key ${params.batch_key} \\
        ${args}

    rm -f merged_for_integration.h5ad

    python ${split_script} \
        --input integrated.h5ad \\
        --sample-key sample_id \\
        --output-dir integrated_samples

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
        harmonypy: \$(python -c "import harmonypy; print(harmonypy.__version__)")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p integration_results integrated_samples
    touch integrated.h5ad
    touch integrated_samples/sample1_integrated.h5ad
    touch integration_results/integration_status.json
    touch versions.yml
    """
}
