process KRAKEN {
    tag "$sample_id"
    container "$params.kraken.docker"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.kraken"), emit: kraken_dir
    tuple val(sample_id), path("${sample_id}.kraken.noise.clean.id"), emit: keep_ids
    path("${sample_id}.report.txt"), emit: report

    
    script:
    """
    kraken2 --db /kraken_db/minikraken2_v1_8GB --paired ${reads[0]} ${reads[1]} --output ${sample_id}.kraken --threads 8 --gzip-compressed --report ${sample_id}.report.txt

    awk '\$3 != "9606" && \$3 !~ /^94[0-9]{2}/ {print \$2}' ${sample_id}.kraken | sed 's/^/@/' > ${sample_id}.kraken.noise.clean.id

    """
}

process SEQTK_PRUNE {
  tag "$sample_id"

  input:
    tuple val(sample_id), path(reads), path(keep_ids)
    
  output:
    tuple val(sample_id), path("${sample_id}.R1.clean.fastq.gz"), path("${sample_id}.R2.clean.fastq.gz")

  script:
  """
  seqtk subseq ${reads[0]} ${keep_ids} | gzip > ${sample_id}.R1.clean.fastq.gz
  seqtk subseq ${reads[1]} ${keep_ids} | gzip > ${sample_id}.R2.clean.fastq.gz
  
  """
}