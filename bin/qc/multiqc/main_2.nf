process MULTIQC_2 {

    tag "Generating MultiQC report"
    
    publishDir "${params.qcdir}/1-QC/genomeQC", mode: 'copy'

    input:
    path quast_dir
    path busco_dir

    output:
    path "multiqc_report"

    script:

    """
    multiqc ${quast_dir} ${busco_dir} -o multiqc_report
    """
}