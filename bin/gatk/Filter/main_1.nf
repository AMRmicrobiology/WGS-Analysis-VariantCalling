process FILTER_VARIANTS {
    tag "Filter Variant ${sample_id}"
    
    publishDir "${params.outdir}", mode: 'copy', saveAs: { filename ->
        if (filename.endsWith(".vcf.gz")) "4-finalVCF/VCF/$filename"
        else null
    }

    container "$params.gatk4.docker"

    input:
    tuple val (sample_id), path(vcf)
    tuple val(id_reference), path(reference)

    output:
    path("${sample_id}_filtered_snp_indel.vcf.gz"), emit: vcf_gz
    tuple val (sample_id), path ("${sample_id}_filtered_snp_indel.vcf.gz"), emit : compl_vcf

    script:
    def referenceBase = reference.baseName
    def referenceDict = referenceBase + ".dict"
    def referenceFai = reference + ".fai"

    """
    echo "Indexing reference ${reference}..."
    if [ ! -f ${referenceFai} ]; then
        samtools faidx ${reference}
    fi

    echo "Creating sequence dictionary for ${reference}..."
    if [ ! -f ${referenceDict} ]; then
        gatk CreateSequenceDictionary -R ${reference} -O ${referenceDict}
    fi

    if [ ! -f ${referenceFai} ]; then
        echo "Error: The reference index (.fai) was not created." >&2
        ls -lh ${reference}
        exit 1
    fi

    if [ ! -f ${referenceDict} ]; then
        echo "Error: The reference dictionary (.dict) was not created." >&2
        ls -lh ${reference}
        exit 1
    fi

    if [ ! -f ${vcf}.tbi ]; then
        echo "Indexing VCF file ${vcf}..."
        tabix -p vcf ${vcf}
    fi
    
    # Filtering SNPs
    gatk VariantFiltration \\
        -R ${reference} \\
        -V ${vcf} \\
        --filter-name "LowQualSNP" --filter-expression "QUAL < 50.0 || MQ < 40.0 || DP < 30" \\
        -O ${sample_id}_snps_filtered.vcf.gz

    # Filtering Indels
    gatk VariantFiltration \\
        -R ${reference} \\
        -V ${vcf} \\
        --filter-name "LowQualIndel" --filter-expression "QUAL < 200.0 || MQ < 40.0 || DP < 30" \\
        -O ${sample_id}_indels_filtered.vcf.gz

    # Select only variants that pass the filter (labels with PASS)
    gatk SelectVariants \\
        -R ${reference} \\
        -V ${sample_id}_snps_filtered.vcf.gz \\
        --exclude-filtered \\
        --select-type-to-include SNP \\
        -O ${sample_id}_snps_pass.vcf.gz

    gatk SelectVariants \\
        -R ${reference} \\
        -V ${sample_id}_indels_filtered.vcf.gz \\
        --exclude-filtered \\
        --select-type-to-include INDEL \\
        -O ${sample_id}_indels_pass.vcf.gz

    # Combine SNPs and indels filtered in one file
    bcftools concat -a -O z -o ${sample_id}_filtered_snp_indel.vcf.gz ${sample_id}_snps_pass.vcf.gz ${sample_id}_indels_pass.vcf.gz
    bcftools index ${sample_id}_filtered_snp_indel.vcf.gz
    """
}