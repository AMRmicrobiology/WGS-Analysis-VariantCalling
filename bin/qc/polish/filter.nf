process FILTER_CONTIGS {
  tag   "FILTER | ${sample_id}"

  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      "docker://${params.short_wgs.docker}" :
      params.short_wgs.docker }"
  
  input:
  tuple val(sample_id), path(contigs)

  output:
  tuple val(sample_id), path("${sample_id}.contigs.filtered.fasta")

  script:
  
  """
  # Contar contigs y longitud total antes de filtrar
  pre_counts=\$(grep -c '^>' ${contigs})
  pre_bases=\$(grep -v '^>' ${contigs} | tr -d '\\n' | wc -c)
  
  # Filtrar contigs <200 nt
  seqtk seq -A -L 300 ${contigs} > ${sample_id}.contigs.filtered.fasta
  
  # Contar después
  post_counts=\$(grep -c '^>' ${sample_id}.contigs.filtered.fasta)
  post_bases=\$(grep -v '^>' ${sample_id}.contigs.filtered.fasta | tr -d '\\n' | wc -c)
  
  # Guardar estadísticas
  cat <<-EOF > filter_stats_${sample_id}.txt
  Sample: ${sample_id}
  Contigs before : \$pre_counts
  Bases before   : \$pre_bases
  Contigs after  : \$post_counts
  Bases after    : \$post_bases
  Removed contigs: \$(( pre_counts - post_counts ))
  Removed bases  : \$(( pre_bases - post_bases ))
  EOF
  """
}