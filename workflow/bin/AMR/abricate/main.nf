process AMR {
    tag "ABRICATE PROCESS"

    publishDir "${params.outdir}/AMR", mode: 'copy'

    container "$params.abricate.docker"

    input:
    tuple val(sample_id), path(assembly_file)

    output:
    path("${sample_id}_abricate_report.tsv"), emit: abricate_report

    script:
    """
    abricate ${assembly_file} > ${sample_id}_abricate_report.tsv
    """
}