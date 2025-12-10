process BUILD_INDEX {
    tag "Index reference genome: ${reference_id.simpleName}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "docker://${params.short_wgs.docker}" :
        params.short_wgs.docker }"
    
    publishDir "${params.reference}/personal/index", mode: 'copy'

    input:
    tuple val(sample_id), path(reference_id)

    output:
    tuple val(sample_id), path("index_${sample_id}.*.bt2"), path(reference_id), emit: index_out

    script:
    """
    bowtie2-build ${reference_id} index_${sample_id}
    """
}