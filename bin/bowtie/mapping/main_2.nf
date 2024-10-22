process ASSEMBLE_GENOME_MAPPING {
    tag "Mapping_Assemble_with_ref. ${sample_id}"
    
    publishDir "${params.outdir}/3-prunning", mode: 'copy', saveAs: { filename ->
        filename.endsWith(".bam") || filename.endsWith(".bai") ? "Pruning_report/$filename" : null
    }

    input:
    tuple val(sample_id), path(asembly)
    path(reference_id)
    
    output:
    tuple val(sample_id), path("${sample_id}.sam"), 
    path("${sample_id}.bam"), path("${sample_id}.bam.bai")

    script:
    """
    bowtie2 -x ${params.index_genome_personal} -f ${asembly} -S ${sample_id}.sam
    samtools sort -o ${sample_id}.bam ${sample_id}.sam
    samtools index ${sample_id}.bam
    """
}