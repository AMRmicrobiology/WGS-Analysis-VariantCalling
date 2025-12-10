process ARIBA {
    tag "ARIBA for ${sample_id} with ${scheme}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    "docker://${params.ariba.docker}" :
    params.ariba.docker }"
    
    publishDir "${params.outdir}/4-MLST/ARIBA", mode: 'copy'

    input:
    tuple val(sample_id), path(trimmed_reads), val(organism), val(scheme)

    output:
    tuple val(sample_id), val(organism), val(scheme), path("outdirectresults/${sample_id}_results_${scheme.replaceAll(/[^a-zA-Z0-9]/, '_')}")

    script:
    """
    # Crear un nombre seguro para el esquema (sin espacios ni caracteres especiales)
    safe_scheme=\$(echo "${scheme}" | sed 's/ /_/g; s/#/_/g')

    # Crear el directorio base `out` si no existe
    mkdir -p out

    # Descargar la base de datos del esquema con ariba pubmlstget
    ariba pubmlstget --verbose "${scheme}" "out/\${safe_scheme}" || { echo "Error: Fallo en la descarga de ${scheme}"; exit 1; }

    # Crear el subdirectorio para los resultados específicos del esquema y la muestra
    mkdir -p "outdirectresults/${sample_id}_results_\${safe_scheme}"

    # Ejecutar ARIBA usando el esquema descargado y las lecturas de la muestra
    ariba run --force "out/\${safe_scheme}/ref_db" \
        "${trimmed_reads[0]}" "${trimmed_reads[1]}" \
        "outdirectresults/${sample_id}_results_\${safe_scheme}" || { echo "Error: Fallo en la ejecución de ARIBA para ${scheme}"; exit 1; }
    """
}