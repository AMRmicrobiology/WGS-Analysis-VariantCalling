process MULTIQC_2 {

    tag "Generating MultiQC report"
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "docker://${params.short_wgs.docker}" :
        params.short_wgs.docker }"
    
    publishDir "${params.outdir}/1-QC/genomeQC", mode: 'copy'

    input:
    path (quast_folder)
    path (busco_folder)

    output:
    path "multiqc_report"

    script:

    """
    multiqc ./ -o multiqc_report

    """
}