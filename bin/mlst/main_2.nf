process MLST {
    tag "MLST-annotation process ${sample_id}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "docker://${params.short_wgs.docker}" :
        params.short_wgs.docker }"
    
    publishDir "${params.outdir}/4-MLST/mlst", mode: 'copy'

    input:
    tuple val(sample_id), path(contigs)

    output:
    path("*.mlst.tab"), emit: 'tab'
    path("*.mlst.json"), emit: 'json'

    script:
    """
    mlst --threads ${task.cpus} --json ${contigs.baseName}.mlst.json ${contigs} > ${contigs.baseName}.mlst.tab
    """
}