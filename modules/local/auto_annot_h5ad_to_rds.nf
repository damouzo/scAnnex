process AUTO_ANNOT_H5AD_TO_RDS {
    tag "global_auto_annot_h5ad_to_rds"
    label 'process_medium'

    conda "\"r-seurat>=5.0\" \"r-base>=4.3\" \"r-reticulate\" \"python\" \"anndata\""
    container "oras://community.wave.seqera.io/library/anndata_python_r-reticulate_r-seurat:7fafb1254a6f3b88"

    input:
    path(h5ad)

    output:
    path "auto_annot_input.rds"               , emit: rds
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    Rscript ${projectDir}/bin/convert_h5ad_to_rds.R \
        --input ${h5ad} \
        --output auto_annot_input.rds \
        --verbose \
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n1 | sed 's/R version //;s/ .*//')
        Seurat: \$(Rscript -e "cat(as.character(packageVersion('Seurat')))" )
    END_VERSIONS
    """

    stub:
    """
    touch auto_annot_input.rds
    touch versions.yml
    """
}
