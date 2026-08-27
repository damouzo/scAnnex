process PSEUDOBULK_BUILD {
    tag "global_pseudobulk_build"
    label 'process_medium'

    // Build pseudo-bulk count matrices from the auto-annotated global object.
    conda "\"scanpy>=1.9\" \"scipy\" \"pandas\""
    container "oras://community.wave.seqera.io/library/scanpy_scipy:af35be00f10024f0"

    input:
    path(h5ad)

    output:
    path "pseudobulk_counts.tsv"        , emit: counts
    path "pseudobulk_metadata.csv"      , emit: metadata
    path "pseudobulk_build_summary.json", emit: build_summary
    tuple path("pseudobulk_counts.tsv"), path("pseudobulk_metadata.csv"), emit: inputs
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    python ${projectDir}/bin/pseudobulk_build.py \
        --input ${h5ad} \
        --counts-out pseudobulk_counts.tsv \
        --metadata-out pseudobulk_metadata.csv \
        --summary-out pseudobulk_build_summary.json \
        --min-cells ${params.pseudobulk_min_cells} \
        --groupby ${params.pseudobulk_groupby} \
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
        anndata: \$(python -c "import anndata; print(anndata.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch pseudobulk_counts.tsv
    touch pseudobulk_metadata.csv
    cat <<-EOF > pseudobulk_build_summary.json
    {
      "success": true,
      "message": "stub"
    }
    EOF
    touch versions.yml
    """
}

process PSEUDOBULK_GLOBAL_BUILD {
    tag "global_pseudobulk_build_all"
    label 'process_medium'

    // Global pseudo-bulk: collapse cell types into a single 'ALL' group per sample.
    // Complements (does not replace) the per-cell-type build. Publish to
    // pseudobulk_dge/global so the two do not collide.
    conda "\"scanpy>=1.9\" \"scipy\" \"pandas\""
    container "oras://community.wave.seqera.io/library/scanpy_scipy:af35be00f10024f0"

    input:
    path(h5ad)

    output:
    path "pseudobulk_global_counts.tsv"        , emit: counts
    path "pseudobulk_global_metadata.csv"      , emit: metadata
    path "pseudobulk_global_build_summary.json", emit: build_summary
    tuple path("pseudobulk_global_counts.tsv"), path("pseudobulk_global_metadata.csv"), emit: inputs
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when
    params.run_pseudobulk_dge

    script:
    def args = task.ext.args ?: ''
    """
    python3 ${projectDir}/bin/pseudobulk_build.py \
        --input ${h5ad} \
        --counts-out pseudobulk_global_counts.tsv \
        --metadata-out pseudobulk_global_metadata.csv \
        --summary-out pseudobulk_global_build_summary.json \
        --min-cells ${params.pseudobulk_min_cells} \
        --groupby ${params.pseudobulk_groupby} \
        --collapse-celltypes \
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
        anndata: \$(python -c "import anndata; print(anndata.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch pseudobulk_global_counts.tsv
    touch pseudobulk_global_metadata.csv
    cat <<-EOF > pseudobulk_global_build_summary.json
    {
      "success": true,
      "message": "stub"
    }
    EOF
    touch versions.yml
    """
}

process PSEUDOBULK_DGE {
    tag "global_pseudobulk_dge"
    label 'process_medium'

    // Differential expression on pseudo-bulk aggregates with DESeq2.
    // Published to GHCR by .github/workflows/build-containers.yml; conda fallback.
    conda "\"bioconductor-deseq2\" \"bioconductor-apeglm\" \"r-optparse\" \"r-jsonlite\" \"r-data.table\" \"r-ggplot2\""
    container "oras://ghcr.io/damouzo/scannex/pseudobulk:latest"

    input:
    path(counts)
    path(metadata)

    output:
    path "pseudobulk_dge/"                      , emit: results_dir
    path "pseudobulk_dge/*_results.csv"         , emit: tables
    path "pseudobulk_dge/*_volcano.{png,pdf}"   , emit: plots, optional: true
    path "pseudobulk_dge/pseudobulk_summary.csv", emit: summary
    path "pseudobulk_dge/pseudobulk_status.json", emit: status_json
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when
    params.run_pseudobulk_dge

    script:
    def args = task.ext.args ?: ''
    def control_group = params.pseudobulk_control_group ?: params.dge_reference
    def shrink = params.pseudobulk_shrink ? "--shrink ${params.pseudobulk_shrink}" : '--shrink apeglm'
    """
    Rscript ${projectDir}/bin/pseudobulk_dge.R \
        --counts ${counts} \
        --metadata ${metadata} \
        --outdir pseudobulk_dge \
        --condition-col ${params.pseudobulk_design} \
        --groupby ${params.pseudobulk_groupby} \
        --sample-col sample_id \
        --control-group ${control_group} \
        --padj-cutoff ${params.pseudobulk_padj_cutoff} \
        ${shrink} \
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n1 | sed 's/R version //;s/ .*//')
        DESeq2: \$(Rscript -e "cat(as.character(packageVersion('DESeq2')))")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p pseudobulk_dge
    touch pseudobulk_dge/contrast_vs_control__celltype_results.csv
    touch pseudobulk_dge/contrast_vs_control__celltype_volcano.png
    touch pseudobulk_dge/pseudobulk_summary.csv
    cat <<-EOF > pseudobulk_dge/pseudobulk_status.json
    {
      "success": true,
      "message": "stub"
    }
    EOF
    touch versions.yml
    """
}

process PSEUDOBULK_GLOBAL_DGE {
    tag "global_pseudobulk_dge_all"
    label 'process_medium'

    // DESeq2 on the global (cell-type-collapsed) pseudo-bulk. Same R script.
    // The R script is pointed at the work directory root so its flat outputs land
    // side by side here and publish cleanly to pseudobulk_dge/global/ without a
    // redundant nested directory.
    conda "\"bioconductor-deseq2\" \"bioconductor-apeglm\" \"r-optparse\" \"r-jsonlite\" \"r-data.table\" \"r-ggplot2\""
    container "oras://ghcr.io/damouzo/scannex/pseudobulk:latest"

    input:
    path(counts)
    path(metadata)

    output:
    path "*_results.csv"                              , emit: tables
    path "*_volcano.{png,pdf}"                        , emit: plots, optional: true
    path "pseudobulk_summary.csv"                     , emit: summary
    path "pseudobulk_status.json"                     , emit: status_json
    path "versions.yml"                               , emit: versions

    when:
    task.ext.when == null || task.ext.when
    params.run_pseudobulk_dge

    script:
    def args = task.ext.args ?: ''
    def control_group = params.pseudobulk_control_group ?: params.dge_reference
    def shrink = params.pseudobulk_shrink ? "--shrink ${params.pseudobulk_shrink}" : '--shrink apeglm'
    """
    Rscript ${projectDir}/bin/pseudobulk_dge.R \
        --counts ${counts} \
        --metadata ${metadata} \
        --outdir . \
        --condition-col ${params.pseudobulk_design} \
        --groupby ${params.pseudobulk_groupby} \
        --sample-col sample_id \
        --control-group ${control_group} \
        --padj-cutoff ${params.pseudobulk_padj_cutoff} \
        ${shrink} \
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version | head -n1 | sed 's/R version //;s/ .*//')
        DESeq2: \$(Rscript -e "cat(as.character(packageVersion('DESeq2')))")
    END_VERSIONS
    """

    stub:
    """
    touch contrast_vs_control__all_results.csv
    touch contrast_vs_control__all_volcano.png
    touch pseudobulk_summary.csv
    cat <<-EOF > pseudobulk_status.json
    {
      "success": true,
      "message": "stub"
    }
    EOF
    touch versions.yml
    """
}