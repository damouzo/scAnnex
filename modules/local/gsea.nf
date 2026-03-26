process GSEA {
    tag "${comparison}"
    label 'process_medium'

    conda "\"bioconductor-clusterprofiler\" \"bioconductor-enrichplot\" \"bioconductor-dose\" \"bioconductor-reactomepa\" \"bioconductor-pathview\" \"bioconductor-keggrest\" \"bioconductor-org.hs.eg.db\" \"bioconductor-org.mm.eg.db\" \"bioconductor-go.db\" \"r-ggplot2\" \"r-data.table\" \"r-dplyr\" \"r-stringr\" \"r-jsonlite\""
    container "docker://damouzo/scannex:gsea-1.0"

    input:
    tuple val(comparison), path(dge_csv)
    val organism
    path gsea_script

    output:
    tuple val(comparison), path("${comparison}_gsea"), emit: results
    path "${comparison}_gsea/*.csv", emit: tables, optional: true
    path "${comparison}_gsea/*.pdf", emit: plots, optional: true
    path "${comparison}_gsea/pathview/*.png", emit: pathview, optional: true
    path "${comparison}_gsea/*.rds", emit: rds_files, optional: true
    path "${comparison}_gsea/*.json", emit: metadata, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def top_pathways = params.gsea_top_pathways ?: 10
    def pathview_top = params.gsea_pathview_top_n ?: 10
    def dashboard_max = params.gsea_dashboard_max_pathways ?: 50
    """
    mkdir -p ${comparison}_gsea/pathview

    Rscript ${gsea_script} \
        --dge-file ${dge_csv} \
        --organism ${organism} \
        --comparison ${comparison} \
        --outdir ${comparison}_gsea \
        --top-pathways ${top_pathways} \
        --pathview-top-n ${pathview_top} \
        --dashboard-max-pathways ${dashboard_max}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(Rscript -e "cat(as.character(getRversion()))")
        clusterProfiler: \$(Rscript -e "cat(as.character(packageVersion('clusterProfiler')))")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${comparison}_gsea/pathview
    touch ${comparison}_gsea/GSEA_GO_BP.csv
    touch ${comparison}_gsea/GSEA_KEGG.csv
    touch ${comparison}_gsea/GSEA_REACTOME.csv
    touch ${comparison}_gsea/gsea_dashboard_data.rds
    touch ${comparison}_gsea/gsea_summary.json
    touch versions.yml
    """
}
