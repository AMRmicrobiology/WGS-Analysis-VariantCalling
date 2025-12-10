process PERSONAL_GENOME_MAPPING {
    tag "Mapping assembly with reference for ${sample_id}"
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "docker://${params.short_wgs.docker}" :
        params.short_wgs.docker }"

    publishDir "${params.outdir}/3-prunning", mode: 'copy', saveAs: { filename ->
        filename.endsWith(".bam") || filename.endsWith(".bai") ? "Pruning_report/$filename" : null
    }

    input:
    tuple val(sample_id), path(reads)
    path(index_mapping_bwa)
    
    output:
    tuple val(sample_id), path("${sample_id}.sam"), 
    path("${sample_id}.bam"), path("${sample_id}.bam.bai"),
    path("${sample_id}_samtools_flagstat.txt")


    script:
    def readGroup = "@RG\\tID:$sample_id\\tSM:$sample_id\\tPL:ILLUMINA"
    
    """
    bwa mem -t ${params.max_threads} -M -R '${readGroup}' ${params.index_mapping_bwa} ${reads[0]} ${reads[1]} > ${sample_id}.sam
    samtools sort < ${sample_id}.sam > ${sample_id}.bam
    samtools index ${sample_id}.bam
    samtools flagstat ${sample_id}.bam > ${sample_id}_samtools_flagstat.txt
    """
}