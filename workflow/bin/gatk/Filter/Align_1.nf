process ALIGN {
    tag "Align Variant ${sample_id}"
    
    container "$params.gatk4.docker"

    input:
    tuple val (sample_id), path(vcf)
    path(reference)

    output:
    tuple val (sample_id), path("${sample_id}_aligned.vcf.gz")

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

    gatk LeftAlignAndTrimVariants \
        -R ${reference} \
        -V ${vcf} \
        -O ${sample_id}_aligned.vcf.gz \
        --split-multi-allelics
    """
}