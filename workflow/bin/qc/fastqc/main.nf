//Quality analisis 
process FASTQC_QUALITY {
    tag "FASTQC"

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