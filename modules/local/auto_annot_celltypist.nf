process AUTO_ANNOT_CELLTYPIST {
    tag "global_auto_annot_celltypist"
    label 'process_medium'

    conda "\"celltypist>=1.6\" \"scanpy>=1.9\""
    container "oras://community.wave.seqera.io/library/celltypist_scanpy:20c2e982b26fecc1"

    input:
    path(h5ad)

    output:
    path "celltypist_annotations.csv"         , emit: annotations
    path "celltypist_status.json"             , emit: status_json
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def models = params.celltypist_models ?: 'Immune_All_Low.pkl'
    def majority_voting = params.celltypist_majority_voting ? '--majority-voting' : ''
    def continue_on_error = params.auto_annot_continue_on_error ? '--continue-on-error' : ''

    """
    python ${projectDir}/bin/auto_annot_celltypist.py \
        --input ${h5ad} \
        --output celltypist_annotations.csv \
        --status celltypist_status.json \
        --models "${models}" \
        ${majority_voting} \
        ${continue_on_error} \
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        celltypist: \$(python -c "import celltypist; print(celltypist.__version__)")
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
    END_VERSIONS
    """

    stub:
    """
    echo "cell_id,auto_annot_celltypist_immune_all_low_pkl,auto_annot_celltypist_immune_all_low_pkl_score" > celltypist_annotations.csv
    echo "cell1,T cells,0.95" >> celltypist_annotations.csv
    cat <<-EOF > celltypist_status.json
    {
      "tool": "celltypist",
      "success": true,
      "message": "stub"
    }
    EOF
    touch versions.yml
    """
}
