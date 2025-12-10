process MARKDUPLICATE {
    tag "Piccar Markduplicate $sample_id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "docker://${params.short_wgs.docker}" :
        params.short_wgs.docker }"

    input:
    tuple val (sample_id), path(bam)

    output:
    tuple val(sample_id), path("${sample_id}.dedup.bam"), emit: dedup_bam
    path("${sample_id}.dedup.metrics.txt")

    script:
    """
    picard MarkDuplicates \
        I=${bam} \
        O=${sample_id}.dedup.bam \
        METRICS_FILE=${sample_id}.dedup.metrics.txt \
        ASSUME_SORT_ORDER=coordinate \
        REMOVE_DUPLICATES=true
    
    samtools index ${sample_id}.dedup.bam
    """
}