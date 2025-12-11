process BAKTA_SET_DB {
    tag "bakta_DB_light"
    
    container {
        workflow.containerEngine == 'singularity' ?
            "docker://${params.bakta_db.docker}" :
            params.bakta_db.docker
    }

    output:
    path("db-light/db-light"), emit: db_bakta_dir

    script:
    """
    export HOME=\$PWD
    DB_ROOT="${params.bakta_db_dir}/db-light/db-light"
    DB_MARKER="version.json"   # fichero que sabemos que está en una DB completa

    if [ ! -f "\${DB_ROOT}/\${DB_MARKER}" ]; then
        echo "[BAKTA_SET_DB] Descargando BD en ${params.bakta_db_dir}/db-light ..."

        set +e
        bakta_db download --type light --output ${params.bakta_db_dir}/db-light
        EXIT=\$?
        set -e

        # Si bakta_db falló pero la DB parece completa (version.json existe), asumimos fallo solo de AMRFinderPlus
        if [ \$EXIT -ne 0 ]; then
            if [ -f "\${DB_ROOT}/\${DB_MARKER}" ]; then
                >&2 echo "[BAKTA_SET_DB] WARNING: bakta_db salió con código \$EXIT; asumo fallo de AMRFinderPlus pero la BD de Bakta está OK."
            else
                >&2 echo "[BAKTA_SET_DB] ERROR: bakta_db falló (código \$EXIT) y no hay \${DB_MARKER}; aborto."
                exit \$EXIT
            fi
        fi
    else
        echo "[BAKTA_SET_DB] Usando BD existente en \${DB_ROOT}"
    fi


    rm -rf db-light
    ln -s ${params.bakta_db_dir}/db-light db-light


    if [ ! -f "db-light/db-light/\${DB_MARKER}" ]; then
        >&2 echo "[BAKTA_SET_DB] ERROR: No encuentro \${DB_MARKER} en db-light/db-light"
        exit 1
    fi
    """
}