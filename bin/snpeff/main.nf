process SNPEFF {
    tag "DB_COMPILATION AND ANNOTATIONS"
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "docker://${params.snpeff.docker}" :
        params.snpeff.docker }"

    publishDir "${params.outdir}/Variant_annotations", mode: 'copy'
    

    input:
    tuple val(new_id), path(variants_vcf), val(id_reference), path(assembly_file), path(gff3_file) ,path(protein_fasta), path(cds_fasta)
    val genome_name_db
    
    output:
    path "annotated_${new_id}_variants.vcf", emit: annotated_vcf

    script:
    """
    
    mkdir -p /opt/conda/envs/snpeff_env/share/snpeff-5.4.0a-0/data/${genome_name_db}

    cp ${assembly_file} /opt/conda/envs/snpeff_env/share/snpeff-5.4.0a-0/data/${genome_name_db}/sequences.fa
    cp ${gff3_file} /opt/conda/envs/snpeff_env/share/snpeff-5.4.0a-0/data/${genome_name_db}/genes.gff
    cp ${protein_fasta} /opt/conda/envs/snpeff_env/share/snpeff-5.4.0a-0/data/${genome_name_db}/protein.fa
    cp ${cds_fasta} /opt/conda/envs/snpeff_env/share/snpeff-5.4.0a-0/data/${genome_name_db}/cds.fa

    echo "# Base de datos para ${genome_name_db}" >> /opt/conda/envs/snpeff_env/share/snpeff-5.4.0a-0/snpEff.config
    echo "${genome_name_db}.genome : ${genome_name_db}" >> /opt/conda/envs/snpeff_env/share/snpeff-5.4.0a-0/snpEff.config

    snpEff build -gff3 -v ${genome_name_db}

    snpEff ann -noLog -noStats -no-upstream -no-downstream -no-utr -c /opt/conda/envs/snpeff_env/share/snpeff-5.4.0a-0/snpEff.config -o vcf ${genome_name_db} ${variants_vcf} > annotated_${new_id}_variants.vcf
    
    """
}