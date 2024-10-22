process AMR_2 {
    tag "AMRFinder PROCESS"

    publishDir "${params.outdir}/AMR/AMRFinder", mode: 'copy'

    container "$params.amrfinderplus.docker"

    input:
    tuple val(sample_id), path(assembly_file)

    output:
    path("${sample_id}_amrfinder_report.tsv"), emit: amrfinder_report

    script:
    """
    amrfinder -n ${assembly_file} -o ${sample_id}_amrfinder_report.tsv
    """
}