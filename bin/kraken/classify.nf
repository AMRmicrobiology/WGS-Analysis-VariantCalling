process KRAKEN2_CLASSIFY {
    tag "${sample_id}"

    // Usa tu imagen con el entrypoint y kraken2 (o la que ya estés usando)
    container "$params.kraken2_new.docker"
    cpus { params.cpus ?: 8 }
    memory { params.memory ?: '16 GB' }
    time '24h'

    // Monta el caché de DB preparado por PREPARE_KRAKEN_DB
    // (en nextflow.config puedes sobreescribir containerOptions si quieres)
    containerOptions "-v ${params.kraken_db_dir}:/kraken2-db"

    /*
     * INPUTS:
     *  - db_ready: sentinel para asegurar que la DB está lista
     *  - tuple(sample_id, [R1, R2?])  ← igual que usas hoy
     */
    input:
    tuple val(sample_id), val(reads)

    /*
     * OUTPUTS:
     *  - raw clasif:  ${sample_id}.kraken
     *  - report:      ${sample_id}.report.txt
     *  - keep_ids:    lista de IDs a conservar (no humanos, ni 94xx)  ← igual a tu awk
     */
    output:
    path "${sample_id}.kraken",        emit: kraken_raw
    path "${sample_id}.report.txt",    emit: report
    path "${sample_id}.keep.ids",      emit: keep_ids

    /*
     * LÓGICA:
     *  - Detecta PE/SE según tamaño de 'reads'
     *  - Añade --gzip-compressed si los ficheros terminan en .gz
     *  - Usa DB en /kraken2-db (preparada por PREPARE_KRAKEN_DB)
     *  - Replica tu awk de filtrado (taxid != 9606 y no 94xx)
     */
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
