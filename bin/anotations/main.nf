process AGAT {
    tag "Merging annotations with AGAT for ${sample_id}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "docker://${params.agat.docker}" :
        params.agat.docker }"
    
    publishDir "${params.outdir}/2-assemble/annotations/data_base", mode: 'copy'

    input:
    path prokka_file
    path bakta_file
    tuple val(sample_id), path(assembly_file)

    output:
    path "fixed_combined_${sample_id}.gff3", emit: combine_gff3
    path "statistics_report_${sample_id}.txt", emit: statistics_report
    path "cds_${sample_id}.fa", emit: cds_fasta
    path "protein_${sample_id}.fa", emit: protein_fasta

    script:
    """
    # Convertir GFF de Prokka a formato GFF3 válido
    agat_convert_sp_gxf2gxf.pl --gff ${prokka_file} --output prokka_${sample_id}.gff3

    # Extraer los contig IDs reales desde el FASTA original
    grep "^>" ${assembly_file} | sed 's/^>//' > bakta_contigs.txt
    grep "^##sequence-region" prokka_${sample_id}.gff3 | cut -d ' ' -f2 > prokka_contigs.txt

    # Validar tamaños iguales
    wc -l prokka_contigs.txt
    wc -l bakta_contigs.txt

    # Generar mapeo (suponiendo orden idéntico de contigs entre ambas anotaciones)
    paste prokka_contigs.txt bakta_contigs.txt > contig_name_map.tsv

    # Renombrar contigs en Prokka
    agat_sq_rename_seqid.pl --gff prokka_${sample_id}.gff3 --tsv contig_name_map.tsv --output prokka_${sample_id}_renamed.gff3

    # Fusionar anotaciones de Prokka y Bakta
    agat_sp_merge_annotations.pl --gff prokka_${sample_id}_renamed.gff3 --gff ${bakta_file} --out combined_${sample_id}.gff3

    # Corregir fases de codón y validar estructura del GFF
    agat_sp_fix_cds_phases.pl --gff combined_${sample_id}.gff3 --fasta ${assembly_file} --output fixed_combined_${sample_id}.gff3

    # Filtrar: conservar solo la isoforma más larga por gen
    agat_sp_keep_longest_isoform.pl --gff fixed_combined_${sample_id}.gff3 --output longest_${sample_id}.gff3

    # Filtrar: eliminar genes sin codificación completa
    agat_sp_filter_incomplete_gene_coding_models.pl --gff longest_${sample_id}.gff3 --fasta ${assembly_file} --output filtered_${sample_id}.gff3

    # Validar consistencia con el FASTA (seqid)
    grep "^>" ${assembly_file} | sed 's/^>//' | sort > ids_fasta.txt
    awk '\$0 ~ /^#/ {next} \$0 ~ /^>/ {exit} {print \$1}' filtered_${sample_id}.gff3 | sort | uniq > ids_gff.txt
    comm -23 ids_gff.txt ids_fasta.txt > mismatched_ids.txt || true

    # Si hay diferencias, renombrar seqid con mapeo invertido
    if [ -s mismatched_ids.txt ]; then
        echo "Hay contigs con nombres incompatibles, ajustando..."
        cut -f2,1 contig_name_map.tsv > map_for_final.tsv
        agat_sq_rename_seqid.pl --gff filtered_${sample_id}.gff3 --tsv map_for_final.tsv --output final_${sample_id}.gff3
    else
        cp filtered_${sample_id}.gff3 final_${sample_id}.gff3
    fi

    # Extraer CDS y proteínas
    gffread final_${sample_id}.gff3 -g ${assembly_file} -x cds_${sample_id}.fa -y protein_${sample_id}.fa

    # Estadísticas de anotación
    agat_sp_statistics.pl --gff final_${sample_id}.gff3 --output statistics_report_${sample_id}.txt
    """
}