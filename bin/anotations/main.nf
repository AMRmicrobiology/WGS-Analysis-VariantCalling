process AGT {
    tag "MERGE ANNOTATIONS"

    publishDir "${params.outdir}/annotations/data_base", mode: 'copy'

    container "$params.agat.docker"
    
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
    agat_convert_sp_gxf2gxf.pl --gff ${prokka_file} --output prokka_${sample_id}.gff3

    agat_sp_merge_annotations.pl --gff prokka_${sample_id}.gff3 --gff ${bakta_file} --out combined_${sample_id}.gff3

    agat_sp_fix_cds_phases.pl --gff combined_${sample_id}.gff3 --fasta ${assembly_file} --output fixed_combined_${sample_id}.gff3

    gffread fixed_combined_${sample_id}.gff3 -g ${assembly_file} -x cds_${sample_id}.fa -y protein_${sample_id}.fa

    agat_sp_statistics.pl --gff fixed_combined_${sample_id}.gff3 --output statistics_report_${sample_id}.txt
    """
}