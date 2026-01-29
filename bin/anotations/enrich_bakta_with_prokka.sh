#!/bin/bash
###############################################################################
# Script: enrich_bakta_with_prokka.sh
# Descripción: Enriquece anotación de BAKTA con información funcional de PROKKA
# Versión: 2.2 - Bacterias | Coordenadas y atributos de Bakta preservados
# Uso: ./enrich_bakta_with_prokka.sh --bakta BAKTA.gff3 --prokka PROKKA.gff \
#                                    --output enriched.gff3 [--min-overlap 0.8]
###############################################################################

set -euo pipefail

# Variables
BAKTA_GFF=""
PROKKA_GFF=""
OUTPUT_GFF=""
MIN_OVERLAP=0.8     # fracción mínima de solapamiento (en Bakta) para enriquecer
VERBOSE=false

log_info() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO: $1"
}

log_error() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERROR: $1" >&2
}

# Parseo de argumentos
while [[ $# -gt 0 ]]; do
    case "$1" in
        --bakta)        BAKTA_GFF="$2"; shift 2 ;;
        --prokka)       PROKKA_GFF="$2"; shift 2 ;;
        --output)       OUTPUT_GFF="$2"; shift 2 ;;
        --min-overlap)  MIN_OVERLAP="$2"; shift 2 ;;
        --verbose)      VERBOSE=true; shift ;;
        -h|--help)
            echo "USO: $0 --bakta FILE --prokka FILE --output FILE [--min-overlap 0.8] [--verbose]"
            exit 0
            ;;
        *)
            log_error "Opción desconocida: $1"
            exit 1
            ;;
    esac
done

# Comprobaciones
if [[ -z "${BAKTA_GFF}" || -z "${PROKKA_GFF}" || -z "${OUTPUT_GFF}" ]]; then
    log_error "Faltan argumentos requeridos (--bakta, --prokka, --output)"
    exit 1
fi

if [[ ! -f "${BAKTA_GFF}" ]]; then
    log_error "Archivo BAKTA no encontrado: ${BAKTA_GFF}"
    exit 2
fi

if [[ ! -f "${PROKKA_GFF}" ]]; then
    log_error "Archivo PROKKA no encontrado: ${PROKKA_GFF}"
    exit 3
fi

log_info "=========================================="
log_info "ENRIQUECIMIENTO BACTERIANO"
log_info "=========================================="
log_info "Base (Bakta):      ${BAKTA_GFF}"
log_info "Fuente (Prokka):   ${PROKKA_GFF}"
log_info "Salida:            ${OUTPUT_GFF}"
log_info "Solapamiento min.: ${MIN_OVERLAP}"
log_info ""

###############################################################################
# PASO 1: Limpiar GFF (cortar en ##FASTA)
###############################################################################

BAKTA_CLEAN="bakta_clean.tmp"
PROKKA_CLEAN="prokka_clean.tmp"

awk '/^##FASTA/ {exit} {print}' "${BAKTA_GFF}" > "${BAKTA_CLEAN}"
awk '/^##FASTA/ {exit} {print}' "${PROKKA_GFF}" > "${PROKKA_CLEAN}"

BAKTA_CDS_COUNT=$(grep -v "^#" "${BAKTA_CLEAN}" | awk '$3=="CDS"' | wc -l)
PROKKA_CDS_COUNT=$(grep -v "^#" "${PROKKA_CLEAN}" | awk '$3=="CDS"' | wc -l)

log_info "CDS en Bakta:  ${BAKTA_CDS_COUNT}"
log_info "CDS en Prokka: ${PROKKA_CDS_COUNT}"
log_info ""

###############################################################################
# PASO 2: Enriquecimiento (AWK con FS=TAB y preservando atributos)
###############################################################################

awk -v min_overlap="${MIN_OVERLAP}" -v verbose="${VERBOSE}" '
BEGIN {
    FS = "\t";    # Separador de ENTRADA: TAB (muy importante)
    OFS = "\t";   # Separador de SALIDA: TAB
    enriched = 0;
    total_cds = 0;
}

########################################################################
# 1) Cargar CDS de PROKKA en memoria (primer fichero: PROKKA_CLEAN)
########################################################################
NR==FNR {
    if ($0 ~ /^#/ || $3 != "CDS") next;

    contig = $1;
    start  = $4;
    end    = $5;
    strand = $7;
    attrs  = $9;

    # Extraer atributos relevantes de Prokka
    product  = "";
    gene     = "";
    ec       = "";
    infer    = "";

    split(attrs, a, ";");
    for (i in a) {
        if (a[i] ~ /^product=/)   product = substr(a[i], 9);
        if (a[i] ~ /^gene=/)      gene    = substr(a[i], 6);
        if (a[i] ~ /^eC_number=/) ec      = substr(a[i], 11);
        if (a[i] ~ /^inference=/) infer   = substr(a[i], 11);
    }

    key = contig ":" start ":" end ":" strand;
    p_product[key]   = product;
    p_gene[key]      = gene;
    p_ec[key]        = ec;
    p_infer[key]     = infer;
    p_list[contig, ++p_count[contig]] = start ":" end ":" strand;

    next;
}

########################################################################
# 2) Procesar BAKTA (segundo fichero: BAKTA_CLEAN)
########################################################################
{
    # Cabeceras, regiones, etc. → se copian tal cual
    if ($0 ~ /^#/)      { print; next }
    if ($3 != "CDS")    { print; next }

    total_cds++;

    b_contig = $1;
    b_start  = $4;
    b_end    = $5;
    b_strand = $7;
    b_attrs  = $9;   # atributos COMPLETOS de Bakta
    b_len    = b_end - b_start + 1;

    best_overlap = 0.0;
    best_product = "";
    best_gene    = "";
    best_ec      = "";
    best_infer   = "";

    # Buscar mejor CDS de Prokka en el mismo contig
    for (i = 1; i <= p_count[b_contig]; i++) {
        split(p_list[b_contig, i], c, ":");
        p_start  = c[1];
        p_end    = c[2];
        p_strand = c[3];

        # Mismo strand
        if (p_strand != b_strand) continue;

        # Solapamiento
        os = (b_start > p_start) ? b_start : p_start;
        oe = (b_end   < p_end)   ? b_end   : p_end;

        if (oe >= os) {
            ol    = oe - os + 1;
            frac  = ol / b_len;
            if (frac > best_overlap) {
                best_overlap = frac;
                key = b_contig ":" p_start ":" p_end ":" p_strand;
                best_product = p_product[key];
                best_gene    = p_gene[key];
                best_ec      = p_ec[key];
                best_infer   = p_infer[key];
            }
        }
    }

    # Si el solapamiento es suficiente, enriquecemos atributos
    if (best_overlap >= min_overlap) {
        enriched++;

        new_attrs = b_attrs;   # arrancamos de los atributos originales de Bakta

        if (best_product != "") new_attrs = new_attrs ";prokka_product="   best_product;
        if (best_gene    != "") new_attrs = new_attrs ";prokka_gene="      best_gene;
        if (best_ec      != "") new_attrs = new_attrs ";prokka_EC_number=" best_ec;
        if (best_infer   != "") new_attrs = new_attrs ";prokka_inference=" best_infer;

        print $1, $2, $3, $4, $5, $6, $7, $8, new_attrs;

        if (verbose == "true") {
            print "[ENRICH] " b_contig ":" b_start "-" b_end " (" best_gene ")" > "/dev/stderr";
        }
    } else {
        # Sin match en Prokka → se deja tal cual
        print;
    }
}

END {
    pct = (total_cds > 0) ? int(enriched * 100 / total_cds) : 0;
    print "" > "/dev/stderr";
    print "========================================" > "/dev/stderr";
    print "RESUMEN DE ENRIQUECIMIENTO" > "/dev/stderr";
    print "========================================" > "/dev/stderr";
    print "Total CDS procesados:   " total_cds  > "/dev/stderr";
    print "CDS enriquecidos:       " enriched   > "/dev/stderr";
    print "Porcentaje enriquecido: " pct "%"    > "/dev/stderr";
    print "========================================" > "/dev/stderr";
}
' "${PROKKA_CLEAN}" "${BAKTA_CLEAN}" > "${OUTPUT_GFF}"

###############################################################################
# PASO 3: Validación final
###############################################################################

log_info "PASO 3: Validando salida..."

if [[ ! -s "${OUTPUT_GFF}" ]]; then
    log_error "No se generó archivo de salida"
    rm -f "${BAKTA_CLEAN}" "${PROKKA_CLEAN}"
    exit 4
fi

FINAL_CDS_COUNT=$(grep -v "^#" "${OUTPUT_GFF}" | awk '$3=="CDS"' | wc -l)
ENRICHED_COUNT=$(grep -v "^#" "${OUTPUT_GFF}" | grep -c "prokka_product=" || echo 0)

log_info "CDS en salida:          ${FINAL_CDS_COUNT}"
log_info "CDS enriquecidos:       ${ENRICHED_COUNT}"

if [[ ${FINAL_CDS_COUNT} -ne ${BAKTA_CDS_COUNT} ]]; then
    log_error "Número de CDS cambió (esperado: ${BAKTA_CDS_COUNT}, obtenido: ${FINAL_CDS_COUNT})"
    rm -f "${BAKTA_CLEAN}" "${PROKKA_CLEAN}"
    exit 5
fi

log_info "✓ Coordenadas de Bakta preservadas"
log_info "✓ Atributos originales de Bakta preservados"
log_info "✓ Anotación funcional enriquecida con Prokka"
log_info "✓ Compatible con SNPeff"

rm -f "${BAKTA_CLEAN}" "${PROKKA_CLEAN}"

log_info "=========================================="
log_info "ENRIQUECIMIENTO COMPLETADO"
log_info "Salida: ${OUTPUT_GFF}"
log_info "=========================================="

exit 0