process AUTO_ANNOT_SCTYPE {
    tag "global_auto_annot_sctype"
    label 'process_medium'

    conda "\"r-base>=4.3\" \"r-optparse\" \"r-jsonlite\" \"r-seurat\" \"r-dplyr\""
    container "oras://community.wave.seqera.io/library/r-azimuth_r-seurat:de4d206e2e153ec1"

    input:
    path(rds)

    output:
    path "sctype_annotations.csv"             , emit: annotations
    path "sctype_status.json"                 , emit: status_json
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def markers_file = params.sctype_markers_file ?: "${projectDir}/assets/sctype_markers_mock.csv"
    def continue_on_error = params.auto_annot_continue_on_error ? '--continue-on-error' : ''
    """
    Rscript ${projectDir}/bin/auto_annot_sctype.R \
        --input ${rds} \
        --output sctype_annotations.csv \
        --status sctype_status.json \
        --markers-file ${markers_file} \
        ${continue_on_error} \
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n1 | sed 's/R version //;s/ .*//')
        Seurat: \$(Rscript -e "cat(as.character(packageVersion('Seurat')))" )
    END_VERSIONS
    """

    stub:
    """
    echo "cell_id,auto_annot_sctype,auto_annot_sctype_score,auto_annot_sctype_cluster" > sctype_annotations.csv
    echo "cell1,T_cells,0.70,cluster_0" >> sctype_annotations.csv
    cat <<-EOF > sctype_status.json
    {
      "tool": "sctype",
      "success": true,
      "message": "stub"
    }
    EOF
    touch versions.yml
    """
}
