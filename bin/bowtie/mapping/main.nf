process PERSONAL_GENOME_MAPPING {
    tag "Mapping personal ref vs ${sample_id}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "docker://${params.short_wgs.docker}" :
        params.short_wgs.docker }"

    cpus 16

    publishDir "${params.outdir}/3-prunning", mode: 'copy',
    saveAs: { filename ->
        filename.endsWith(".bam") || filename.endsWith(".bai") ? "Pruning_report/$filename" : null
    }

    input:
    tuple val(sample_id), path(reads), val(ref_id), path(index_files), path(reference_fasta)

    output:
    tuple val(sample_id), 
          path("${sample_id}.sam"), 
          path("${sample_id}.sorted.bam"),
          path("${sample_id}.sorted.bam.bai"),
          path("${sample_id}_bowtie2_metrics.txt"),
          path("${sample_id}_samtools_flagstat.txt"), emit: all_outputs

    script:
    def index_prefix = "index_${ref_id}"

    """
    # Mapping con Bowtie2
    bowtie2 -x ${index_prefix} -1 ${reads[0]} -2 ${reads[1]} -S ${sample_id}.sam \
      --very-sensitive \
      -D 20 -R 3 -L 22 -i S,1,0.50 \
      --mp 6,2 \
      --np 1 \
      --rdg 5,3 \
      --rfg 5,3 \
      -p ${task.cpus} \
      --no-unal \
      --met-file ${sample_id}_bowtie2_metrics.txt

    # Procesamiento del SAM
    samtools sort -o ${sample_id}.sorted.bam ${sample_id}.sam
    samtools index ${sample_id}.sorted.bam
    samtools flagstat ${sample_id}.sorted.bam > ${sample_id}_samtools_flagstat.txt
    """
}
