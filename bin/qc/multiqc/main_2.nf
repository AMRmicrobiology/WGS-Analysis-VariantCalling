process MULTIQC_2 {

    tag "Generating MultiQC report"
    
    publishDir "${params.outdir}/1-QC/genomeQC", mode: 'copy'

    input:
    tuple val (sample_id), path (quast_dir)
    tuple val (sample_id), path (busco_dir)

    output:
    path "multiqc_report"

    script:

    """
    multiqc ${params.quast_dir} ${params.busco_dir} -o multiqc_report
    """
}