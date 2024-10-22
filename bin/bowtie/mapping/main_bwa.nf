process ASSEMBLE_GENOME_MAPPING {
    tag "Mapping_Assemble_with_ref. ${sample_id}"
    
    publishDir "${params.outdir}/3-prunning", mode: 'copy', saveAs: { filename ->
        filename.endsWith(".bam") || filename.endsWith(".bai") ? "Pruning_report/$filename" : null
    }

    input:
    tuple val(sample_id), path(asembly)
    path(params.personal_ref)
    
    output:
    tuple val(sample_id), path("${sample_id}.sam"), 
    path("${sample_id}.bam"), path("${sample_id}.bam.bai")

    script:
    def readGroup = "@RG\\tID:$sample_id\\tSM:$sample_id\\tPL:ILLUMINA"
    
    """
    bwa mem -t ${params.max_threads} -M -R '$readGroup' ${params.personal_ref} ${asembly} > ${sample_id}.sam
    samtools sort < ${sample_id}.sam > ${sample_id}.bam
    samtools index ${sample_id}.bam

    """
}