//Quality analisis 
process FASTQC_QUALITY {
    tag "FASTQC"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "docker://${params.short_wgs.docker}" :
        params.short_wgs.docker }"

    input:
    path (reads)

    output:
    path ("*.html"), emit:qc_html
    path ("*.zip"), emit: qc_zip

    script:
    """
    fastqc ${reads}
    """
}