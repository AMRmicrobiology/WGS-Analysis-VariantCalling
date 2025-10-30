process FILTER_VARIANTS {
    tag "Filter Variant ${sample_id}"
    
    publishDir "${params.outdir}", mode: 'copy', saveAs: { filename ->
        if (filename.endsWith(".vcf.gz")) "4-VCF/filter_VCF/$filename"
        else null
    }

    container "$params.gatk4.docker"

    input:
    tuple val (sample_id), path(vcf), val(id_reference), path(reference)

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

    # Filtering SNPs with Balance of Alleles
    echo "Filtering SNPs with Quality and Balance of Alleles..."
    gatk VariantFiltration \\
        -R ${reference} \\
        -V ${vcf} \\
        --filter-name "LowQualSNP" \\
        --filter-expression "QUAL < 100.0 || MQ < 40.0 || DP < 50 || QD < 2.0 || FS > 60.0 || SOR > 3.0 || (vc.getGenotype(0).getAD() == null || (vc.getGenotype(0).getAD().1 + 1.0) / (vc.getGenotype(0).getAD().0 + vc.getGenotype(0).getAD().1 + 1.0) < 0.95)" \\
        -O ${sample_id}_snps_filtered.vcf.gz

    # Filtering Indels with Homopolymer Regions
    echo "Filtering Indels with Quality and Homopolymer Regions..."
    gatk VariantFiltration \\
        -R ${reference} \\
        -V ${vcf} \\
        --filter-name "LowQualIndel" \\
        --filter-expression "QUAL < 200.0 || MQ < 40.0 || DP < 30 || QD < 2.0 || FS > 200.0 || SOR > 10.0 || HRun > 6" \\
        -O ${sample_id}_indels_filtered.vcf.gz

    # Select only variants that pass the filter (labels with PASS)
    echo "Selecting passing SNPs..."
    gatk SelectVariants \\
        -R ${reference} \\
        -V ${sample_id}_snps_filtered.vcf.gz \\
        --exclude-filtered \\
        --select-type-to-include SNP \\
        -O ${sample_id}_snps_pass.vcf.gz

    echo "Selecting passing Indels..."
    gatk SelectVariants \\
        -R ${reference} \\
        -V ${sample_id}_indels_filtered.vcf.gz \\
        --exclude-filtered \\
        --select-type-to-include INDEL \\
        -O ${sample_id}_indels_pass.vcf.gz

    # Combine SNPs and indels filtered in one file
    echo "Combining SNPs and Indels..."
    bcftools concat -a -O z -o ${sample_id}_filtered_snp_indel.vcf.gz ${sample_id}_snps_pass.vcf.gz ${sample_id}_indels_pass.vcf.gz
    bcftools index ${sample_id}_filtered_snp_indel.vcf.gz
    """
}