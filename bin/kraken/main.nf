process KRAKEN {
    tag "$sample_id"
    container "$params.kraken.docker"

    cpus   { params.kraken_cpus }
    memory { params.kraken_mem  }
    time '24h'

    input:
    tuple val(sample_id), path(reads), path db_dir

    output:
    tuple val(sample_id), path("${sample_id}.kraken"), emit: kraken_dir
    tuple val(sample_id), path("${sample_id}.kraken.noise.clean.id"), emit: keep_ids
    path("${sample_id}.report.txt"), emit: report

    script:
    
    """
    kraken2 \
    --db "${db_dir}" \
    --paired "${reads[0]}" "${reads[1]}" \
    --threads ${task.cpus} \
    --gzip-compressed \
    --memory-mapping \
    ${ params.kraken2_extra_args ?: '' } \
    ${ params.kraken_confidence ? "--confidence ${params.kraken_confidence}" : "" } \
    --use-names \
    --report "${sample_id}.report.txt" \
    > "${sample_id}.kraken"

    awk '\$3 != "9606" && \$3 !~ /^94[0-9]{2}/ {print \$2}' ${sample_id}.kraken > ${sample_id}.kraken.noise.clean.id
    """
}

process SEQTK_PRUNE {
  tag "$sample_id"

  input:
    tuple val(sample_id), path(reads), path(keep_ids)
    
  output:
    tuple val(sample_id), path("${sample_id}.R{1,2}.clean.fastq.gz"), emit: pruned_reads

  script:
  """
  seqtk subseq ${reads[0]} ${keep_ids} | gzip > ${sample_id}.R1.clean.fastq.gz
  seqtk subseq ${reads[1]} ${keep_ids} | gzip > ${sample_id}.R2.clean.fastq.gz
  
  """
}