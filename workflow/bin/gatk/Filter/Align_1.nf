process ALIGN {
    tag "Align Variant ${sample_id}"
    
    container "$params.gatk4.docker"

    input:
    tuple val (sample_id), path(vcf)
    tuple val(id_reference), path(reference)

    output:
    tuple val (sample_id), path("${sample_id}_aligned.vcf.gz")

    script:
    def referenceDict = reference.toString().replaceAll('\\.(fna|fa|fasta)$', '.dict')
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
        exit 1
    fi

    if [ ! -f ${referenceDict} ]; then
        echo "Error: The reference dictionary (.dict) was not created." >&2
        exit 1
    fi

    if [ ! -f ${vcf}.tbi ]; then
        echo "Indexing VCF file ${vcf}..."
        gatk IndexFeatureFile -I ${vcf}
    fi

    echo "Running LeftAlignAndTrimVariants"
    gatk LeftAlignAndTrimVariants \
        -R ${reference} \
        -V ${vcf} \
        -O ${sample_id}_aligned.vcf.gz \
        --split-multi-allelics
    """
}