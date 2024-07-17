process PERSONAL_GENOME_MAPPING {
    tag "Mapping_personal_ref ${sample_id}"
    
    publishDir "${params.outdir}/3-prunning", mode: 'copy',
    saveAs: { filename ->
        filename.endsWith(".bam") || filename.endsWith(".bai") ? "Pruning_report/$filename" : null
    }

    input:
    tuple val(sample_id), path(reads)
    path (reference_id)
    
    output:
    tuple val(sample_id), path("${sample_id}.sam"), 
    path("${sample_id}.bam"), path("${sample_id}.bam.bai"),
    path("${sample_id}_bowtie2_metrics.txt"), path("${sample_id}_samtools_flagstat.txt")
    
    
    script:
    """
    bowtie2 -x ${params.index_genome_personal} -1 ${reads[0]} -2 ${reads[1]} -S ${sample_id}.sam \
    --very-sensitive \
    --met-file ${sample_id}_bowtie2_metrics.txt
    samtools sort < ${sample_id}.sam > ${sample_id}.bam
    samtools index ${sample_id}.bam
    samtools flagstat ${sample_id}.bam > ${sample_id}_samtools_flagstat.txt
    """
}