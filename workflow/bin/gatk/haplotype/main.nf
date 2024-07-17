process HAPLOTYPECALLER {
    tag "Haplotype ${sample_id}"

    container "$params.gatk4.docker"

    input:
    tuple val(sample_id), path(bam)
    path(reference)

    output:
    path("${sample_id}.g.vcf.gz")

    
    script:
    def referenceDict = reference.toString().replaceAll('\\.fna$', '.dict')
    
    """
    # Ensure the reference is indexed
    samtools faidx ${reference}
    
    # Create the dictionary if it does not exist
    if [ ! -f ${referenceDict} ]; then
        gatk CreateSequenceDictionary -R ${reference} -O ${referenceDict}
    fi
    
    # Check if BAM is indexed, if not, index it
    if [ ! -f ${bam}.bai ]; then
        samtools index ${bam}
    fi
    
    gatk --java-options "-Xmx4g" HaplotypeCaller \
        --sample-ploidy 1 \
        -R ${reference} \
        -I ${bam} \
        -O ${sample_id}.g.vcf.gz \
        --standard-min-confidence-threshold-for-calling 30 \
        --minimum-mapping-quality 20 \
        -ERC GVCF
    """
}
