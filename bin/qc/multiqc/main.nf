process MULTIQC {

    tag "Generating MultiQC report"

    publishDir "${params.outdir}/1-QC/fastqQC", mode: 'copy'

    input:
    path fastqc_first
    path fastqc_after

    output:
    path "multiqc_report"

    script:
    """
    echo "FastQC files: ${fastqc_first} ${fastqc_after}"

    multiqc ${fastqc_first} ${fastqc_after} -o multiqc_report
    """
}