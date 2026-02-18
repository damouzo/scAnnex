process AUTO_ANNOT_CELLTYPIST {
    tag "$meta.id"
    label 'process_medium'

    conda "\"celltypist>=1.6\" \"scanpy>=1.9\""
    container "oras://community.wave.seqera.io/library/celltypist_scanpy:20c2e982b26fecc1"

    input:
    tuple val(meta), path(h5ad)

    output:
    tuple val(meta), path("*_annotated.h5ad"), emit: h5ad
    path "*_celltypist.csv"                   , emit: annotations
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    // SLC: Build CellTypist parameters
    def model = params.celltypist_model ?: 'Immune_All_Low.pkl'
    def majority_voting = params.celltypist_majority_voting ? '--majority-voting' : ''
    
    """
    # Run CellTypist annotation and save to H5AD
    auto_annot_celltypist.py \\
        --input ${h5ad} \\
        --output ${prefix}_celltypist.csv \\
        --output-h5ad ${prefix}_annotated.h5ad \\
        --model ${model} \\
        ${majority_voting} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        celltypist: \$(python -c "import celltypist; print(celltypist.__version__)")
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_annotated.h5ad
    echo "cell_id,label,score,tool" > ${prefix}_celltypist.csv
    echo "cell1,T cells,0.95,celltypist" >> ${prefix}_celltypist.csv
    touch versions.yml
    """
}
