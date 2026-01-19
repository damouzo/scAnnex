process H5AD_TO_RDS {
    tag "$meta.id"
    label 'process_medium'

    container "docker.io/satijalab/seurat:5.0.0"

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("*.rds"), emit: rds
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    convert_h5ad_to_rds.R \\
        --input ${h5ad} \\
        --output ${prefix}.rds \\
        --verbose \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n1 | sed 's/R version //;s/ .*//')
        Seurat: \$(Rscript -e "cat(as.character(packageVersion('Seurat')))")
        SeuratDisk: \$(Rscript -e "cat(as.character(packageVersion('SeuratDisk')))")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.rds
    touch versions.yml
    """
}
