process SNPEFF {

    tag "DB_COMPILATION AND ANNOTATIONS"

    publishDir "${params.outdir}/annotations", mode: 'copy'
    
    container "$params.snpeff.docker"

    input:
    path gff3_file
    tuple val(sample_id), path(assembly_file)
    val genome_name_db
    path protein_fasta
    path cds_fasta
    tuple val(new_id), path(variants_vcf)

    output:
    path "annotated_${new_id}_variants.vcf", emit: annotated_vcf

    script:
    """
    # Crear el directorio para el nuevo genoma dentro del directorio correcto de SnpEff
    mkdir -p /opt/conda/envs/snpeff_env/share/snpeff-5.2-1/data/${genome_name_db}

    # Copiar el archivo FASTA y GFF al nuevo directorio, renombrando el GFF3 a GFF
    cp ${assembly_file} /opt/conda/envs/snpeff_env/share/snpeff-5.2-1/data/${genome_name_db}/sequences.fa
    cp ${gff3_file} /opt/conda/envs/snpeff_env/share/snpeff-5.2-1/data/${genome_name_db}/genes.gff
    cp ${protein_fasta} /opt/conda/envs/snpeff_env/share/snpeff-5.2-1/data/${genome_name_db}/protein.fa
    cp ${cds_fasta} /opt/conda/envs/snpeff_env/share/snpeff-5.2-1/data/${genome_name_db}/cds.fa

    # Modificar snpeff.config para agregar el nuevo genoma al final del archivo
    echo "# Base de datos para ${genome_name_db}" >> /opt/conda/envs/snpeff_env/share/snpeff-5.2-1/snpEff.config
    echo "${genome_name_db}.genome : ${genome_name_db}" >> /opt/conda/envs/snpeff_env/share/snpeff-5.2-1/snpEff.config

    # Compilar la base de datos para SnpEff
    snpEff build -gff3 -v ${genome_name_db}


    # Ejecutar SnpEff para anotar el archivo VCF descomprimido
    snpEff ann -noLog -noStats -no-upstream -no-downstream -no-utr -c /opt/conda/envs/snpeff_env/share/snpeff-5.2-1/snpEff.config -o vcf ${genome_name_db} ${variants_vcf} > annotated_${new_id}_variants.vcf
    """
}