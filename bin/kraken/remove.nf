process KRAKENTOOLS_EXCLUDE {
  tag { sample_id }

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "docker://${params.kraken_tools.docker}" :
        params.kraken_tools.docker }"

  input:
    tuple val(sample_id), path(reads), path(kraken_out)
  
  output:
    tuple val(sample_id), path("${sample_id}.R{1,2}.clean.fastq.gz"), emit: pruned_reads
  
  script:
  
  """
  extract_kraken_reads.py \
    -k ${kraken_out} \
    -s ${reads[0]} -s2 ${reads[1]} \
    -o ${sample_id}.R1.clean.fq -o2 ${sample_id}.R2.clean.fq \
    --exclude --taxid 9443 --include-children
  gzip ${sample_id}.R1.clean.fq ${sample_id}.R2.clean.fq
  """
}