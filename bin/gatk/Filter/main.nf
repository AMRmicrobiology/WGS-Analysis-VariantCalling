process FILTER_VARIANTS {
    tag "Filter Variant ${type}"
    
    publishDir "${params.outdir}", mode: 'copy', saveAs: { filename ->
        if (filename.endsWith(".vcf.gz")) "4-finalVCF/VCF/$filename"
        else null
    }

    container "$params.gatk4.docker"

    input:
    path(vcf)
    path(reference)
    val type

    output:
    path("${type}_filtered_snp_indel.vcf.gz")

    script:
    def referenceDict = reference.toString().replaceAll('\\.(fna|fa)$', '.dict')
    
    """
    # Ensure the reference is indexed
    samtools faidx ${reference}
    
    # Create the dictionary if it does not exist
    if [ ! -f ${referenceDict} ]; then
        gatk CreateSequenceDictionary -R ${reference} -O ${referenceDict}
    fi

    # Indexar el archivo VCF de entrada si no está indexado
    if [ ! -f ${vcf}.tbi ]; then
        gatk IndexFeatureFile -I ${vcf}
    fi

    # Filtrar SNPs
    gatk VariantFiltration \\
        -R ${reference} \\
        -V ${vcf} \\
        --filter-name "LowQualSNP" --filter-expression "QUAL < 50.0 || MQ < 25.0 || DP < 30" \\
        -O ${type}_snps_filtered.vcf.gz

    # Filtrar Indels
    gatk VariantFiltration \\
        -R ${reference} \\
        -V ${vcf} \\
        --filter-name "LowQualIndel" --filter-expression "QUAL < 200.0 || MQ < 25.0 || DP < 30" \\
        -O ${type}_indels_filtered.vcf.gz

    # Seleccionar solo variantes que pasen los filtros (etiquetadas como PASS)
    gatk SelectVariants \\
        -R ${reference} \\
        -V ${type}_snps_filtered.vcf.gz \\
        --exclude-filtered \\
        --select-type-to-include SNP \\
        -O ${type}_snps_pass.vcf.gz

    gatk SelectVariants \\
        -R ${reference} \\
        -V ${type}_indels_filtered.vcf.gz \\
        --exclude-filtered \\
        --select-type-to-include INDEL \\
        -O ${type}_indels_pass.vcf.gz

    # Combinar SNPs e indels filtrados en un solo archivo
    bcftools concat -a -O z -o ${type}_filtered_snp_indel.vcf.gz ${type}_snps_pass.vcf.gz ${type}_indels_pass.vcf.gz
    bcftools index ${type}_filtered_snp_indel.vcf.gz
    """
}