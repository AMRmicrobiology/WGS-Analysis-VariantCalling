process GENOTYPE {
    tag "genotype ${sample_id}"
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "docker://${params.gatk4.docker}" :
        params.gatk4.docker }"

    publishDir "${params.outdir}", mode: 'copy', saveAs: { filename ->
        if (filename.endsWith(".vcf.gz")) "4-VCF/genotype_VCF/$filename"
        else null
    }
    
    input:
    tuple val(sample_id), path(vcf), val(id_reference), path(reference)

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
        --sample-ploidy 1 \
        --max-alternate-alleles 6 \
        --allow-old-rms-mapping-quality-annotation-data false \
        --annotate-with-num-discovered-alleles true

    # Index the output VCF file
    tabix -f -p vcf final_${sample_id}.vcf.gz

    echo "Generating VCF statistics with bcftools..."
    bcftools stats final_${sample_id}.vcf.gz > final_${sample_id}.vcf.stats
    """
}