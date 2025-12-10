process ALIGMENT_PILON {
  tag "ALIGMENT | ${sample_id}"
  cpus 8
  memory '32 GB'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "docker://${params.short_wgs.docker}" :
        params.short_wgs.docker }"

  input:
  tuple val(sample_id), path(filtered_fasta), path(reads_clear)

  output:
  tuple val(sample_id), path("${sample_id}.bwa.aln.sorted.bam"), emit: aln_bam
  path("${sample_id}.bwa.aln.sorted.bam.bai"), emit: aln_bai
  path("${sample_id}.bwa.aln.sorted.bam.stats.txt"), emit: aln_stats

  script:
  """
  # 1) Indexa el ensamblaje
  bwa index ${filtered_fasta}
  # 2) Alinea lecturas y ordena
  bwa mem -t ${task.cpus} ${filtered_fasta} ${reads_clear[0]} ${reads_clear[1]} \
    | samtools sort -O BAM -o ${sample_id}.bwa.aln.sorted.bam
  samtools index ${sample_id}.bwa.aln.sorted.bam

  # 3) Estadísticas del alineamiento
  samtools flagstat ${sample_id}.bwa.aln.sorted.bam > ${sample_id}.bwa.aln.sorted.bam.stats.txt
  samtools idxstats ${sample_id}.bwa.aln.sorted.bam >> ${sample_id}.bwa.aln.sorted.bam.stats.txt
  
  """
}

process PILON_POLISH {
  tag   "PILON | ${sample_id}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "docker://${params.pilon.docker}" :
        params.pilon.docker }"

  cpus 8
  memory '32 GB'

  publishDir "${params.outdir}/2-Assembly/", mode: 'copy', patterns: ['*.pilon.changes.txt', '${sample_id}.fasta']

  input:
  tuple val(sample_id), path(filtered_fasta), path(index_bam)

  output:
  tuple val(sample_id), path("${sample_id}.fasta"), emit: pilon_fa
  path "${sample_id}.pilon.changes.txt", emit: info_pilon_changes

  script:
  """

  java -Xmx${task.memory.toGiga()}g \
       -jar /pilon/pilon.jar \
       --genome ${filtered_fasta} \
       --frags ${index_bam} \
       --output pilon \
       --threads ${task.cpus} \
       --changes \
       --vcf \
       --tracks

  # Move and save stats results
  mv pilon.fasta ${sample_id}.fasta
  mv pilon.changes ${sample_id}.pilon.changes.txt
  mv pilon.vcf ${sample_id}.pilon.variants.vcf

  """
}