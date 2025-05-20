process SNPEFF {

    tag "DB_COMPILATION AND ANNOTATIONS"
    publishDir "${params.outdir}/Variant_annotations", mode: 'copy'
    container "$params.snpeff.docker"
    cpus 2

    input:
    path gff3_file
    tuple val(id_reference), path(assembly_file)
    val genome_name_db
    tuple val(new_id), path(variants_vcf)


    output:
    path "annotated_${new_id}_variants.vcf", emit: annotated_vcf

    script:
    """
    # Variables de entorno
    SNPEFF_HOME=/opt/conda/envs/snpeff_env/share/snpeff-5.2-1
    DATA_DIR=\\\$SNPEFF_HOME/data

    # 1) Prepara la DB
    mkdir -p \$DATA_DIR/${genome_name_db}
    cp ${assembly_file} \$DATA_DIR/${genome_name_db}/sequences.fa
    cp ${gff3_file}     \$DATA_DIR/${genome_name_db}/genes.gff

    echo "${genome_name_db}.genome : ${genome_name_db}" \
      >> \$SNPEFF_HOME/snpEff.config

    # 2) Construye la base de datos
    snpEff build -gff3 -c \$SNPEFF_HOME/snpEff.config -noCheckCds -noCheckProtein ${genome_name_db}

    snpEff ann -c \$SNPEFF_HOME/snpEff.config -noLog -noStats -no-upstream -no-downstream -no-utr -v ${genome_name_db} ${variants_vcf} > annotated_${new_id}_variants.vcf
    """
}