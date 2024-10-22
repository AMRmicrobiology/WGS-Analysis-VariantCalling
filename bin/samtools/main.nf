process MERGE_BAM {
    tag "merge BAMs"


    input:
    path(bam)

    output:
    tuple val ("merge"), path ("merged_replicas.bam"), emit: merge_data
    path ("merged_replicas.bam.bai")

    script:
    """
    samtools merge merged_replicas.bam ${bam.join(' ')}
    samtools index merged_replicas.bam
    """
}