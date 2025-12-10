process MULTIQC {
    tag "Generating MultiQC report"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "docker://${params.short_wgs.docker}" :
        params.short_wgs.docker }"
    
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