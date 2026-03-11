process AUTO_ANNOT_SUMMARIZE {
    tag "global_auto_annot_summarize"
    label 'process_medium'

    conda "\"scanpy>=1.9\" \"pandas\""
    container "oras://community.wave.seqera.io/library/scanpy_scipy:af35be00f10024f0"

    input:
    path(base_h5ad)
    path(celltypist_annotations)
    path(celltypist_status)
    path(sctype_annotations)
    path(sctype_status)
    path(azimuth_annotations)
    path(azimuth_status)
    path(singler_annotations)
    path(singler_status)

    output:
    path "auto_annotated_global.h5ad"         , emit: h5ad
    path "auto_annot_summary.json"            , emit: summary
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    python ${projectDir}/bin/auto_annot_merge.py \
        --input ${base_h5ad} \
        --output auto_annotated_global.h5ad \
        --summary auto_annot_summary.json \
        --celltypist ${celltypist_annotations} \
        --celltypist-status ${celltypist_status} \
        --sctype ${sctype_annotations} \
        --sctype-status ${sctype_status} \
        --azimuth ${azimuth_annotations} \
        --azimuth-status ${azimuth_status} \
        --singler ${singler_annotations} \
        --singler-status ${singler_status} \
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch auto_annotated_global.h5ad
    cat <<-EOF > auto_annot_summary.json
    {
      "success": true,
      "message": "stub"
    }
    EOF
    touch versions.yml
    """
}
