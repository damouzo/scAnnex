process AUTO_ANNOT_AZIMUTH {
    tag "global_auto_annot_azimuth"
    label 'process_medium'

    conda "\"r-azimuth=0.5.0\" \"r-seurat=5.4.0\" \"r-optparse\" \"r-jsonlite\""
    container "oras://community.wave.seqera.io/library/r-azimuth_r-seurat:de4d206e2e153ec1"

    input:
    path(rds)

    output:
    path "azimuth_annotations.csv"            , emit: annotations
    path "azimuth_status.json"                , emit: status_json
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def refs = params.azimuth_refs ?: 'Human - PBMC'
    def continue_on_error = params.auto_annot_continue_on_error ? '--continue-on-error' : ''
    """
    Rscript ${projectDir}/bin/auto_annot_azimuth.R \
        --input ${rds} \
        --output azimuth_annotations.csv \
        --status azimuth_status.json \
        --refs "${refs}" \
        ${continue_on_error} \
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n1 | sed 's/R version //;s/ .*//')
        Azimuth: \$(Rscript -e "cat(as.character(packageVersion('Azimuth')))" )
        Seurat: \$(Rscript -e "cat(as.character(packageVersion('Seurat')))" )
    END_VERSIONS
    """

    stub:
    """
    echo "cell_id,auto_annot_azimuth_human_pbmc,auto_annot_azimuth_human_pbmc_score,auto_annot_azimuth_human_pbmc_l1,auto_annot_azimuth_human_pbmc_l2" > azimuth_annotations.csv
    echo "cell1,T cell,0.88,Lymphoid,T" >> azimuth_annotations.csv
    cat <<-EOF > azimuth_status.json
    {
      "tool": "azimuth",
      "success": true,
      "message": "stub"
    }
    EOF
    touch versions.yml
    """
}
