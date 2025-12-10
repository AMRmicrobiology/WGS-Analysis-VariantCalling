process AGT {
    tag "Merging annotations with AGAT for ${sample_code}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "docker://${params.agat.docker}" :
        params.agat.docker }"

    publishDir "${params.outdir}/2-Assembly/3-Annotations/AGT_${sample_code}", mode: 'copy'

    input:
    path prokka_file
    path bakta_file
    tuple val(sample_code), path(assembly_file)

    output:
    path "final_${sample_code}.gff3", emit: final_gff3
    path "statistics_report_${sample_code}.txt", emit: statistics_report
    path "cds_${sample_code}.fa", emit: cds_fasta
    path "protein_${sample_code}.fa", emit: protein_fasta
    path "QC_report_${sample_code}.txt", emit: qc_report

    script:
    """
    set -Eeuo pipefail
    echo "🧬 ===== [AGT] Fusión de anotaciones para ${sample_code} ====="

    ########################################
    # 1. Validar archivos de entrada
    ########################################
    for f in "${prokka_file}" "${bakta_file}" "${assembly_file}"; do
        if [[ ! -s "\$f" ]]; then
            echo "❌ [ERROR] Archivo de entrada faltante o vacío: \$f" >&2
            exit 1
        fi
    done

    echo "[1/12] Limpieza de FASTA (IDs únicos por contig)..."
    awk '/^>/{print ">contig_" ++i; next} {print}' ${assembly_file} > assembly_clean.fasta

    ########################################
    # 2. Convertir GFF de Prokka a GFF3
    ########################################
    echo "[2/12] Convirtiendo anotación de Prokka..."
    agat_convert_sp_gxf2gxf.pl --gff ${prokka_file} --output prokka_${sample_code}.gff3

    ########################################
    # 3. Validar GFFs
    ########################################
    echo "[3/12] Validando archivos GFF..."
    agat_validate_gff.pl --gff prokka_${sample_code}.gff3 > validate_prokka.log || true
    agat_validate_gff.pl --gff ${bakta_file} > validate_bakta.log || true

    ########################################
    # 4. Mapeo 1:1 entre contigs originales y renombrados
    ########################################
    echo "[4/12] Generando mapeo 1:1 de contigs..."
    grep "^>" ${assembly_file} | sed 's/^>//; s/ .*//' > prokka_contigs.txt
    grep "^>" assembly_clean.fasta | sed 's/^>//' > bakta_contigs.txt

    total_prokka=\$(wc -l < prokka_contigs.txt)
    total_bakta=\$(wc -l < bakta_contigs.txt)
    if [[ "\$total_prokka" -ne "\$total_bakta" ]]; then
        echo "⚠️  [WARNING] Diferente número de contigs (\$total_prokka vs \$total_bakta). Se usará mapeo por orden." >&2
    fi

    paste prokka_contigs.txt bakta_contigs.txt > contig_name_map.tsv

    ########################################
    # 5. Renombrar contigs
    ########################################
    echo "[5/12] Renombrando contigs en anotación Prokka..."
    agat_sq_rename_seqid.pl --gff prokka_${sample_code}.gff3 --tsv contig_name_map.tsv --output prokka_${sample_code}_renamed.gff3

    ########################################
    # 6. Fusionar anotaciones Prokka + Bakta
    ########################################
    echo "[6/12] Fusionando anotaciones..."
    agat_sp_merge_annotations.pl --gff prokka_${sample_code}_renamed.gff3 --gff ${bakta_file} --out combined_${sample_code}.gff3

    ########################################
    # 7. Corrección de fases y limpieza
    ########################################
    echo "[7/12] Corrigiendo fases de codón..."
    agat_sp_fix_cds_phases.pl --gff combined_${sample_code}.gff3 --fasta assembly_clean.fasta --output fixed_combined_${sample_code}.gff3

    echo "[8/12] Conservando isoforma más larga..."
    agat_sp_keep_longest_isoform.pl --gff fixed_combined_${sample_code}.gff3 --output longest_${sample_code}.gff3

    echo "[9/12] Filtrando genes incompletos..."
    agat_sp_filter_incomplete_gene_coding_models.pl --gff longest_${sample_code}.gff3 --fasta assembly_clean.fasta --output filtered_${sample_code}.gff3

    ########################################
    # 10. Validar integridad de IDs y coordenadas
    ########################################
    echo "[10/12] Validando consistencia de contig IDs..."
    grep "^>" assembly_clean.fasta | sed 's/^>//' | sort > ids_fasta.txt
    awk '\$0 !~ /^#/ {print \$1}' filtered_${sample_code}.gff3 | sort | uniq > ids_gff.txt
    comm -23 ids_gff.txt ids_fasta.txt > mismatched_ids.txt || true

    if [ -s mismatched_ids.txt ]; then
        echo "⚠️  Contigs incompatibles detectados, ajustando..."
        cut -f2,1 contig_name_map.tsv > map_for_final.tsv
        agat_sq_rename_seqid.pl --gff filtered_${sample_code}.gff3 --tsv map_for_final.tsv --output final_${sample_code}.gff3
    else
        cp filtered_${sample_code}.gff3 final_${sample_code}.gff3
    fi

    ########################################
    # 11. QC bacteriano adicional
    ########################################
    echo "[11/12] Ejecución de controles de calidad..."
    {
        echo "### QC Report for ${sample_code}"
        echo "Genome size: \$(awk '/^[^>]/ {sum+=length(\$0)} END{print sum}' assembly_clean.fasta) bp"
        echo "Number of contigs: \$(grep -c '^>' assembly_clean.fasta)"
        echo "Number of CDS: \$(grep -c 'CDS' final_${sample_code}.gff3)"
        echo "Number of genes: \$(grep -c 'gene' final_${sample_code}.gff3)"
        echo "Number of exons: \$(grep -c 'exon' final_${sample_code}.gff3)"
        echo "Number of RNAs: \$(grep -c 'rRNA\\|tRNA' final_${sample_code}.gff3)"
        echo ""
        echo "# Checking gene structure (bacterial sanity checks)"
        awk '\$3=="mRNA" && \$5<\$4 {print "❌ Gene with inverted coordinates:", \$0}' final_${sample_code}.gff3 | tee -a QC_report_${sample_code}.txt || true
        awk '\$3=="exon" && \$5<\$4 {print "❌ Exon with invalid coordinates:", \$0}' final_${sample_code}.gff3 | tee -a QC_report_${sample_code}.txt || true
        echo "✅ QC check completed."
    } > QC_report_${sample_code}.txt

    ########################################
    # 12. Extracción de CDS y estadísticas
    ########################################
    echo "[12/12] Extrayendo CDS y proteínas..."
    gffread final_${sample_code}.gff3 -g assembly_clean.fasta -x cds_${sample_code}.fa -y protein_${sample_code}.fa

    echo "Generando estadísticas..."
    agat_sp_statistics.pl --gff final_${sample_code}.gff3 --output statistics_report_${sample_code}.txt

    echo "===== [AGT] Finalizado correctamente para ${sample_code} ====="
    """
}