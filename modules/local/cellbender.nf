process CELLBENDER {
    tag "$meta.id"
    label 'process_high'

    conda "\"scanpy>=1.9\" \"cellbender>=0.3.0\" \"pytables\""
    container "oras://community.wave.seqera.io/library/pytables_scanpy_pip_cellbender:6fa233dea96313cf"

    input:
    tuple val(meta), path(raw_matrix), path(filtered_matrix)

    output:
    tuple val(meta), path("*_cellbender.h5ad"), emit: h5ad
    path "cellbender_output/*.h5"              , emit: h5_output
    path "cellbender_output/*.pdf"             , emit: plots
    path "cellbender_output/*.log"             , emit: log
    path "versions.yml"                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    
    // Build CellBender arguments
    def expected_cells = params.cellbender_expected_cells ? "--expected-cells ${params.cellbender_expected_cells}" : ''
    def total_droplets = params.cellbender_total_droplets ? "--total-droplets ${params.cellbender_total_droplets}" : ''
    def fpr = params.cellbender_fpr ? "--fpr ${params.cellbender_fpr}" : '--fpr 0.01'
    def epochs = params.cellbender_epochs ? "--epochs ${params.cellbender_epochs}" : '--epochs 150'
    
    """
    cellbender_remove_background.py \\
        --raw-matrix ${raw_matrix} \\
        --filtered-matrix ${filtered_matrix} \\
        --output ${prefix}_cellbender.h5ad \\
        --output-dir cellbender_output \\
        ${expected_cells} \\
        ${total_droplets} \\
        ${fpr} \\
        ${epochs} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellbender: \$(python -c "import cellbender; print(cellbender.__version__)" 2>/dev/null || echo "0.3.0")
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p cellbender_output
    touch ${prefix}_cellbender.h5ad
    touch cellbender_output/cellbender_output.h5
    touch cellbender_output/cellbender_output.pdf
    touch cellbender_output/cellbender_output.log
    touch versions.yml
    """
}
