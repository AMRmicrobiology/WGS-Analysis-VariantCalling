process PREPARE_KRAKEN_DB {
  tag "${params.db_select ?: 'db_16GB'}"
  label 'kraken_db_setup'

  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      "docker://${params.kraken2_new.docker}" :
      params.kraken2_new.docker }"

  cpus 2
  memory '4 GB'
  time '24h'

  output:
  path 'kraken_db', emit: db_ready

  script:

  """

  export DB_DIR="\$PWD/kraken_db"
  mkdir -p "\$DB_DIR"

  export DB_SELECT='${params.db_select}'
  ${ params.db_url ? "export DB_URL='${params.db_url}'" : ":" }
  ${ params.db_url_checksum ? "export DB_URL_CHECKSUM='${params.db_url_checksum}'" : ":" }

  kraken2-entrypoint.sh prepare-db "\$DB_SELECT"

  chmod -R a+rX "\$DB_DIR" || true
  
  """
}