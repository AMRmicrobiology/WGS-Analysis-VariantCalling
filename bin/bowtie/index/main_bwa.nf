process BUILD_INDEX_1 {
    tag "Indexing ${sample_id}"
    label 'index_process'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "docker://${params.short_wgs.docker}" :
        params.short_wgs.docker }"

    
    publishDir "${params.reference}/personal", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reference_id)

    output:
    path("*")

    script:
    """
	bwa index ${reference_id}
    """
}