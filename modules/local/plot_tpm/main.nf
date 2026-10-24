process PLOT_TPM {
    label 'process_single'

    // built with Seqera Containers (pandas, matplotlib, procps-ng)
    container 'community.wave.seqera.io/library/matplotlib_pandas_procps-ng:4f598c8cdc89f96c'

    input:
    path tpm                    // tpm_matrix.tsv from TPM_TABLE

    output:
    path "*_mqc.png", emit: plots
    tuple val("${task.process}"), val('matplotlib'), eval("python3 -c 'import matplotlib; print(matplotlib.__version__)'"), topic: versions, emit: versions_matplotlib

    script:
    """
    plot_tpm.py ${tpm}
    """
}