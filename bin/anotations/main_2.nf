process ENRICHMENT_ANNOTATION {
    tag "Enriching annotations with AGAT for ${sample_id}"
    label 'annotation'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "docker://jimmlucas/enrichment:v1.0.0" :
        "jimmlucas/enrichment:v1.0.0" }"

    publishDir "${params.outdir}/2-Assembly/3-Annotations/ENRICHED_${sample_id}", mode: 'copy'
    
    input:
    tuple val(sample_id), path(prokka_file), path(bakta_file), path(assembly_file)

    output:
    path "enriched_${sample_id}.gff3", emit: enriched_gff3
    tuple val(sample_id), path("enriched_${sample_id}.gff3"), emit: enriched_gff3_tuple
    path "cds_${sample_id}.fa", emit: cds_fasta
    path "protein_${sample_id}.fa", emit: protein_fasta
    path "enrichment_report_${sample_id}.txt", emit: report

    script:
    """
    set -euo pipefail

    # 1) Inicializar reporte
    echo "========================================" > enrichment_report_${sample_id}.txt
    echo "ENRIQUECIMIENTO DE ANOTACIÓN BACTERIANA" >> enrichment_report_${sample_id}.txt
    echo "========================================" >> enrichment_report_${sample_id}.txt
    
    # 2) Preparar script de enriquecimiento
    cp ${projectDir}/bin/anotations/enrich_bakta_with_prokka.sh .
    chmod +x enrich_bakta_with_prokka.sh
    
    # 3) Ejecutar enriquecimiento (Mantiene coordenadas de Bakta)
    ./enrich_bakta_with_prokka.sh \
        --bakta ${bakta_file} \
        --prokka ${prokka_file} \
        --output enriched_${sample_id}.gff3 \
        --verbose >> enrichment_report_${sample_id}.txt 2>&1
    
    # 4) Limpieza específica para bacterias antes de gffread
    # - Eliminamos líneas con strand '?' (como oriC) que rompen gffread
    # - Eliminamos tabs accidentales en la columna de atributos que causan warnings
    echo "[INFO] Limpiando GFF para extracción de secuencias..." >> enrichment_report_${sample_id}.txt
    
    grep -v \$'\t?\t' enriched_${sample_id}.gff3 | sed 's/\t/ /g9' > enriched_clean.gff3

    # 5) Extraer CDS y Proteínas (FASTA)
    # Usamos -y para proteínas y -x para CDS, esencial para validación de SNPs
    echo "[INFO] Generando archivos FASTA..." >> enrichment_report_${sample_id}.txt
    gffread enriched_clean.gff3 \
        -g ${assembly_file} \
        -x cds_${sample_id}.fa \
        -y protein_${sample_id}.fa

    echo "[INFO] Proceso bacteriano completado." >> enrichment_report_${sample_id}.txt
    """
}