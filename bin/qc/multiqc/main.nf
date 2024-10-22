process MULTIQC {

    tag "Generating MultiQC report"
    
    publishDir "${params.qcdir}", mode: 'copy'

    input:
    path fastqc_first
    path fastqc_after
    tuple val(sample_id), path(quast)

    output:
    path "multiqc_report"
    path "multiqc_report_assemble"

    script:

    """
    multiqc ${fastqc_first} ${fastqc_after} -o multiqc_report

    multiqc ${quast} -o multiqc_report_assemble
    """
}