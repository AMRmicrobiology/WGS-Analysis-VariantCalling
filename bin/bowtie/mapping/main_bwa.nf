process PERSONAL_GENOME_MAPPING {
    tag "Mapping ${sample_id}"

    publishDir "${params.outdir}/3-prunning", mode: 'copy', saveAs: { filename ->
        filename.endsWith(".bam") || filename.endsWith(".bai") ? "Pruning_report/$filename" : null
    }

    input:
    tuple val(sample_id), path(reads), path(fasta_index), path(index_reference_bwa)

    output:
    tuple val(sample_id), 
          path("${sample_id}.sam"), 
          path("${sample_id}.bam"), 
          path("${sample_id}.bam.bai"), 
          path("${sample_id}_samtools_flagstat.txt")

    script:
    def readGroup = "@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA"

    """
    # Mapeo con BWA-MEM
    bwa mem -t ${params.max_threads} -M -R '${readGroup}' ${fasta_index} ${reads[0]} ${reads[1]} > ${sample_id}.sam

    # Conversión a BAM ordenado
    samtools sort -o ${sample_id}.bam ${sample_id}.sam

    # Indexación BAM
    samtools index ${sample_id}.bam

    # Estadísticas de mapeo
    samtools flagstat ${sample_id}.bam > ${sample_id}_samtools_flagstat.txt
    """
}