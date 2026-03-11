process AUTO_ANNOT_SINGLER {
    tag "global_auto_annot_singler"
    label 'process_medium'

    conda "\"r-base>=4.3\" \"bioconductor-singler\" \"bioconductor-celldex\" \"bioconductor-singlecellexperiment\" \"r-seurat\""
    container "oras://community.wave.seqera.io/library/bioconductor-celldex_bioconductor-singlecellexperiment_bioconductor-singler_r-seurat:f688d516296a8a25"

    input:
    path(rds)

    output:
    path "singler_annotations.csv"            , emit: annotations
    path "singler_status.json"                , emit: status_json
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def refs = params.singler_refs ?: 'BlueprintEncodeData'
    def prune = params.singler_prune ? '--prune' : ''
    def continue_on_error = params.auto_annot_continue_on_error ? '--continue-on-error' : ''
    """
    Rscript ${projectDir}/bin/auto_annot_singler.R \
        --input ${rds} \
        --output singler_annotations.csv \
        --status singler_status.json \
        --refs "${refs}" \
        ${prune} \
        ${continue_on_error} \
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n1 | sed 's/R version //;s/ .*//')
        SingleR: \$(Rscript -e "cat(as.character(packageVersion('SingleR')))" )
        celldex: \$(Rscript -e "cat(as.character(packageVersion('celldex')))" )
    END_VERSIONS
    """

    stub:
    """
    echo "cell_id,auto_annot_singler_blueprintencodedata,auto_annot_singler_blueprintencodedata_score,auto_annot_singler_blueprintencodedata_delta_next,auto_annot_singler_blueprintencodedata_pruned" > singler_annotations.csv
    echo "cell1,T_cells,0.80,0.12,T_cells" >> singler_annotations.csv
    cat <<-EOF > singler_status.json
    {
      "tool": "singler",
      "success": true,
      "message": "stub"
    }
    EOF
    touch versions.yml
    """
}
