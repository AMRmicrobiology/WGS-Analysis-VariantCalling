process SNPEFF {
    tag { "SNPEFF_${new_id}" }

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "docker://${params.snpeff.docker}" :
        params.snpeff.docker }"

    publishDir "${params.outdir}/Variant_annotations", mode: 'copy'
    
    input:
    tuple(
        path(gff3_file),
        val(id_reference),
        path(assembly_file),
        val(genome_name_db),
        val(new_id),
        path(variants_vcf)
    )

    output:
    tuple val(new_id), path("annotated_${new_id}_variants.vcf"), emit: annotated_vcf_tuple
    path("annotated_${new_id}_variants.vcf"), emit: annotated_vcf

    script:
    """
    echo ">>> GFF3: ${gff3_file}"
    echo ">>> REF ID: ${id_reference}"
    echo ">>> Ensamblado FASTA: ${assembly_file}"
    echo ">>> DB nombre: ${genome_name_db}"
    echo ">>> Sample ID: ${new_id}"
    echo ">>> VCF input: ${variants_vcf}"

    SNPEFF_HOME=/opt/conda/envs/snpeff_env/share/snpeff-5.2-1
    DB_DIR=\$PWD/snpeff_db/${genome_name_db}
    CONFIG_FILE=\$PWD/snpEff.config

    mkdir -p \$DB_DIR

    # Copiar los archivos necesarios solo si no existen
    [ ! -f "\$DB_DIR/sequences.fa" ] && cp ${assembly_file} \$DB_DIR/sequences.fa
    [ ! -f "\$DB_DIR/genes.gff" ]     && cp ${gff3_file}     \$DB_DIR/genes.gff

    # Configuración de snpEff
    if [ ! -f "\$CONFIG_FILE" ]; then
      echo "${genome_name_db}.genome : ${genome_name_db}" > \$CONFIG_FILE
    fi

    # Construcción condicional de la base de datos
    if [ ! -f "\$DB_DIR/snpEffectPredictor.bin" ]; then
      echo "→ Construyendo DB SNPeff..."
      snpEff build -gff3 -c \$CONFIG_FILE -dataDir \$PWD/snpeff_db -noCheckCds -noCheckProtein ${genome_name_db}
    else
      echo "→ DB SNPeff ya existe. Saltando construcción."
    fi

    # Anotación
    snpEff ann -c \$CONFIG_FILE -dataDir \$PWD/snpeff_db -noLog -noStats -no-upstream -no-downstream -no-utr -v ${genome_name_db} ${variants_vcf} > annotated_${new_id}_variants.vcf
    """
}