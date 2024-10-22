process BUILD_INDEX_1 {
    tag "index"
    label 'index_process'
   
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