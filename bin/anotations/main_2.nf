process AGT {
    tag "Merging annotations with AGAT for ${sample_id}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "docker://${params.agat.docker}" :
        params.agat.docker }"

    publishDir "${params.outdir}/2-assemble/annotations/data_base", mode: 'copy'

    input:
    tuple val(sample_id), path(prokka_file), path(bakta_file)
    tuple val(sample_id), path(assembly_file)

    output:
    path "fixed_combined_${sample_id}.gff3", emit: combine_gff3
    path "statistics_report_${sample_id}.txt", emit: statistics_report
    path "cds_${sample_id}.fa", emit: cds_fasta
    path "protein_${sample_id}.fa", emit: protein_fasta

    script:
    """
    # 1. Convertir el GFF de Prokka a GFF3 válido
    agat_convert_sp_gxf2gxf.pl --gff ${prokka_file} --output prokka_${sample_id}.gff3

    # 2. Extraer nombres de contigs desde Prokka y Bakta para crear el mapeo
    grep -v '^#' prokka_${sample_id}.gff3 | cut -f1 | sort | uniq > prokka_contigs.txt
    grep -v '^#' ${bakta_file} | cut -f1 | sort | uniq > bakta_contigs.txt

    # 3. Generar archivo de mapeo de nombres (asumiendo orden de contigs idéntico)
    paste prokka_contigs.txt bakta_contigs.txt > contig_name_map.tsv

    # 4. Renombrar los contigs del GFF de Prokka para que coincidan con Bakta
    agat_sq_rename_seqid.pl --gff prokka_${sample_id}.gff3 --tsv contig_name_map.tsv --output prokka_${sample_id}_renamed.gff3

    # 5. Fusionar anotaciones de Prokka y Bakta
    agat_sp_merge_annotations.pl --gff prokka_${sample_id}_renamed.gff3 --gff ${bakta_file} --out combined_${sample_id}.gff3

    # 6. Corregir fases de codón y validar estructura del GFF
    agat_sp_fix_cds_phases.pl --gff combined_${sample_id}.gff3 --fasta ${assembly_file} --output fixed_combined_${sample_id}.gff3

    # 7. Extraer CDS y proteínas
    gffread fixed_combined_${sample_id}.gff3 -g ${assembly_file} -x cds_${sample_id}.fa -y protein_${sample_id}.fa

    # 8. Generar informe estadístico
    agat_sp_statistics.pl --gff fixed_combined_${sample_id}.gff3 --output statistics_report_${sample_id}.txt
    """
}