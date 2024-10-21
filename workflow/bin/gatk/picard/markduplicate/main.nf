process MARKDUPLICATE {
    tag "Piccar Markduplicate $sample_id"
    

    input:
    tuple val (sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.dedup.bam"), path("${sample_id}.dedup.metrics.txt")

    script:
    """
    picard MarkDuplicates \
        I=${bam} \
        O=${sample_id}.dedup.bam \
        METRICS_FILE=${sample_id}.dedup.metrics.txt \
        ASSUME_SORTED=True \
        REMOVE_DUPLICATES=true
    """
}