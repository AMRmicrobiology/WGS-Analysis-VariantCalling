process BUILD_INDEX_1 {
    tag "index"
    label 'index_process'
   
    publishDir "${params.reference}/personal/index_bwa", mode: 'copy'
    
    input:
    tuple val(sample_id), path(reference_id)

    output:
    tuple val(sample_id), path(reference_id), emit: fasta
    path("*"), emit: index_files

    script:
    """
	bwa index ${reference_id}
    """
}