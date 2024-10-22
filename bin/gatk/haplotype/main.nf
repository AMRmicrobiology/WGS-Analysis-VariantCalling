process HAPLOTYPECALLER {
    tag "Haplotype ${sample_id}"

    container "$params.gatk4.docker"

    input:
    tuple val(sample_id), path(bam)
    path(reference)

    output:
    path("${sample_id}.g.vcf.gz")

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

    # Verificar que los archivos se hayan creado
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

    echo "Indexing BAM file ${bam}..."
    if [ ! -f ${bam}.bai ]; then
        samtools index ${bam}
    fi

    # Verificar que el archivo BAM esté indexado
    if [ ! -f ${bam}.bai ]; then
        echo "Error: The BAM index (.bai) was not created." >&2
        exit 1
    fi

    echo "Running HaplotypeCaller for ${sample_id}..."
    gatk --java-options "-Xmx4g" HaplotypeCaller \
        --sample-ploidy 1 \
        -R ${reference} \
        -I ${bam} \
        -O ${sample_id}.g.vcf.gz \
        --standard-min-confidence-threshold-for-calling 30 \
        --minimum-mapping-quality 20 \
        -ERC GVCF

    # Verificar que el archivo de salida se haya creado
    if [ ! -f ${sample_id}.g.vcf.gz ]; then
        echo "Error: The output VCF file (.g.vcf.gz) was not created." >&2
        exit 1
    fi
    """
}