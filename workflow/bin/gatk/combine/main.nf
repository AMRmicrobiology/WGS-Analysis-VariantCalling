process COMBINE_GVCFS {
    tag "CombineGVCFs"

    container "$params.gatk4.docker"

    input:
    path gvcf_files
    path reference
    val type

    output:
    path "combined_${type}.vcf.gz"

    script:
    def referenceDict = reference.toString().replaceAll('\\.fna$', '.dict')

    """
    # Ensure the reference is indexed
    if [ ! -f ${reference}.fai ]; then
        echo "Indexing reference fasta"
        samtools faidx ${reference}
    fi

    # Create the dictionary if it does not exist
    if [ ! -f ${referenceDict} ]; then
        echo "Creating reference dictionary"
        gatk CreateSequenceDictionary -R ${reference} -O ${referenceDict}
    fi

    # Ensure all GVCF files are indexed
    for gvcf in ${gvcf_files}; do
        if [ ! -f \${gvcf}.tbi ]; then
            echo "Indexing GVCF file: \${gvcf}"
            gatk IndexFeatureFile -I \${gvcf}
        fi
    done

    # Combine GVCFs
    gatk CombineGVCFs \\
        -R ${reference} \\
        ${gvcf_files.collect { "--variant " + it }.join(' \\\n        ')} \\
        -O combined_${type}.vcf.gz
    """
}