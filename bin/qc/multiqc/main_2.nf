process MULTIQC_2 {

    tag "Generating MultiQC report"
    
    publishDir "${params.outdir}/1-QC/genomeQC", mode: 'copy'

    input:
    tuple val (sample_id), path (quast_folder)
    tuple val (sample_id), path (busco_folder)

    output:
    path "multiqc_report"

    script:

    """
    multiqc ./ -o multiqc_report

    """
}