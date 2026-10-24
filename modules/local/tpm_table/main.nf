process TPM_TABLE {
    label 'process_single'

    // MultiQC image: contains Python 3 and 'ps' (needed by Nextflow); the script only uses the Python standard library
    container 'community.wave.seqera.io/library/multiqc:1.35--c17fb751507e9dfc'

    input:
    path counts                 // all featureCounts tables (collected)

    output:
    path "tpm_matrix.tsv", emit: tpm
    tuple val("${task.process}"), val('python'), eval("python3 --version | sed 's/Python //'"), topic: versions, emit: versions_python

    script:
    """
    counts_to_tpm.py ${counts} > tpm_matrix.tsv
    """
}