//Quality analisis 
process FASTQC_QUALITY {
    tag "FASTQC"

    publishDir "${params.qcdir}", mode: 'copy'

    input:
    path (reads)

    output:
    path ("*.html")

    script:
    """
    fastqc ${reads}
    """
}