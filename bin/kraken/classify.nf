process KRAKEN2_CLASSIFY {
    tag "${sample_id}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "docker://${params.kraken2_new.docker}" :
        params.kraken2_new.docker }"

    cpus { params.cpus ?: 8 }
    memory { params.memory ?: '16 GB' }
    time '24h'

    containerOptions "-v ${params.kraken_db_dir}:/kraken2-db"

    input:
    tuple val(sample_id), val(reads)


    output:
    path "${sample_id}.kraken",        emit: kraken_raw
    path "${sample_id}.report.txt",    emit: report
    path "${sample_id}.keep.ids",      emit: keep_ids


    script:
    def r1      = reads[0]
    def r2      = (reads.size() > 1 ? reads[1] : null)
    def gzFlag  = reads.every{ it.toString().endsWith('.gz') } ? '--gzip-compressed' : ''
    def extra   = params.kraken2_extra_args ?: ''
    def dbArg   = '--db /kraken2-db'     // si tu módulo ya pone --db, puedes quitar esta línea

    if (r2) {
        """
        kraken2 ${dbArg} --threads ${task.cpus} ${gzFlag} ${extra} \
            --paired ${r1} ${r2} \
            --output ${sample_id}.kraken \
            --report ${sample_id}.report.txt

        # Igual a tu awk original: columna 3 = taxid; columna 2 = read-id
        awk '\$3 != "9606" && \$3 !~ /^94[0-9]{2}/ {print \$2}' ${sample_id}.kraken > ${sample_id}.keep.ids
        """
    } else {
        """
        kraken2 ${dbArg} --threads ${task.cpus} ${gzFlag} ${extra} \
            ${r1} \
            --output ${sample_id}.kraken \
            --report ${sample_id}.report.txt

        awk '\$3 != "9606" && \$3 !~ /^94[0-9]{2}/ {print \$2}' ${sample_id}.kraken > ${sample_id}.keep.ids
        """
    }
}
