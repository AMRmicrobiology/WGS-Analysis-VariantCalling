process BUILD_INDEX {
    tag "index"
    label 'index_process'
   
    publishDir "${params.reference}/personal", mode: 'copy'
    
    input:
    path reference_id

    output:
    path("*")

    script:
    """
	bwa index ${reference_id}
    """
}