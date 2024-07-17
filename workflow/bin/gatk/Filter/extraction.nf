process SELECT_ANALYSY {
    tag "selection of ${type}"
    
    publishDir "${params.outdir}", mode: 'copy', saveAs: { filename ->
        if (filename.endsWith(".vcf.gz")) "4-finalVCF/VCF/$filename"
        else null
    }

    container "$params.gatk4.docker"

    input:
    path(vcf)
    path(reference)
    val type

    output:
    path("${type}_wildtype_variants_AB1.vcf.gz")
    path ("${type}_treated_samples_variants.vcf.gz")

    script:
    def referenceDict = reference.toString().replaceAll('\\.fna$', '.dict')
    
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
    
    gatk SelectVariants \
    -R ${reference} \
    -V ${vcf} \
    -sn AB1 \
    -O ${type}_wildtype_variants_AB1.vcf.gz

    gatk SelectVariants \
    -R ${reference} \
    -V ${vcf} \
    -sn AB2 -sn AB3 -sn AB4 -sn AB5 -sn AB6 -sn AB7 -sn AB8 -sn AB9 -sn AB10 -sn AB11 \
    -O ${type}_treated_samples_variants.vcf.gz
    """
}