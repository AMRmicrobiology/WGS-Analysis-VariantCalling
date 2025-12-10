process BUSCO {
    tag "GENOME COMPLETENESS ${sample_id}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "docker://${params.busco.docker}" :
        params.busco.docker }"
    
    publishDir "${params.outdir}/1-QC/genomeQC/BUSCO", mode: "copy" 

    input:
    tuple val(sample_id), path(assemble)

    output:
    tuple val(sample_id), path("${sample_id}_busco")

    script:

    """
    busco -i ${assemble} -m genome -l bacteria -o ${sample_id}_busco
    """
}
