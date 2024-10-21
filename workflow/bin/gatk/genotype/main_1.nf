process GENOTYPE {
    tag "genotype ${sample_id}"
    
    publishDir "${params.outdir}", mode: 'copy', saveAs: { filename ->
        if (filename.endsWith(".vcf.gz")) "4-finalVCF/VCF/$filename"
        else null
    }
    
    container "$params.gatk4.docker"  

    input:
    tuple val(sample_id), path(vcf)
    tuple val(sample_id), path(reference)

    output:
    tuple val(sample_id), path("final_${sample_id}.vcf.gz")

    script:
    def referenceDict = reference.toString().replaceAll('\\.(fna|fa)$', '.dict')
    
    """
    if [ ! -f ${reference}.fai ]; then
        echo "Indexing reference fasta"
        samtools faidx ${reference}
    fi
    if [ ! -f ${referenceDict} ]; then
        echo "Creating reference dictionary"
        gatk CreateSequenceDictionary -R ${reference} -O ${referenceDict}
    fi
    if [ ! -f ${vcf}.tbi ]; then
        echo "Indexing vcf file"
        tabix -p vcf ${vcf}
    fi
    
    echo "Running GenotypeGVCFs"

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