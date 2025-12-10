process KRAKEN {
  tag "$sample_id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "docker://${params.kraken.docker}" :
        params.kraken.docker }"


  cpus   { params.kraken_cpus }
  memory { params.kraken_mem  }
  time '24h'

  input:
  tuple val(sample_id), path(reads), path (db_dir)

  output:
  tuple val(sample_id), path("${sample_id}.kraken"), emit: kraken_dir
  tuple val(sample_id), path("${sample_id}.kraken.noise.clean.id"), emit: keep_ids
  path("${sample_id}.report.txt"), emit: report

  script:
    """
    kraken2 \\
      --db "${db_dir}" \\
      --paired "${reads[0]}" "${reads[1]}" \\
      --threads ${task.cpus} \\
      --gzip-compressed \\
      --memory-mapping \\
      ${ params.kraken2_extra_args ?: '' } \\
      ${ params.kraken_confidence ? "--confidence ${params.kraken_confidence}" : "" } \\
      --report "${sample_id}.report.txt" \\
      > "${sample_id}.kraken"

    # keep_ids excluyendo Primates y ancestros (taxid ID)
    
    cat > ids.awk << 'AWK'
    BEGIN{
      split("9443 9606 9605 9604 9598 9593 9601 9526 9483 314295 40674", a, " ");
      for(i in a) deny[a[i]]=1;
    }
    {
      status=\$1; rid=\$2; tax=\$3;
      # U = Unclassified -> conservar
      if (status=="U") { print rid; next }
      # C = Classified -> conservar solo si NO está en la denylist
      if (!(tax in deny)) { print rid }
    }
    AWK

    awk -f ids.awk ${sample_id}.kraken > ${sample_id}.kraken.noise.clean.id
    """
}


process SEQTK_PRUNE {
  tag "$sample_id"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "docker://${params.short_wgs.docker}" :
        params.short_wgs.docker }"

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