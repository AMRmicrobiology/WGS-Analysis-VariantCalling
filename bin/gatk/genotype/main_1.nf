process GENOTYPE {
    tag "genotype ${sample_id}"
    
    publishDir "${params.outdir}", mode: 'copy', saveAs: { filename ->
        if (filename.endsWith(".vcf.gz")) "4-finalVCF/VCF/$filename"
        else null
    }
    
    container "$params.gatk4.docker"  

    input:
    tuple val(sample_id), path(vcf)
    tuple val(id_reference), path(reference)

    output:
    tuple val(sample_id), path("final_${sample_id}.vcf.gz")

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

    gatk GenotypeGVCFs \
        -R ${reference} \
        -V ${vcf} \
        -O final_${sample_id}.vcf.gz \
        --max-alternate-alleles 6 \
        --allow-old-rms-mapping-quality-annotation-data false

    # Index the output VCF file
    tabix -f -p vcf final_${sample_id}.vcf.gz
    """
}