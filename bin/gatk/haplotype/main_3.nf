process HAPLOTYPECALLER {
    tag "Haplotype ${sample_id}"
    
    publishDir "${params.outdir}", mode: 'copy', saveAs: { filename ->
        if (filename.endsWith(".vcf.")) "8-tryVCF/VCF/$filename"
        else null
    }

    container "$params.gatk4.docker"

    input:
    tuple val (sample_id), path (bam)
    tuple val(barcode), path (reference)

    output:
    tuple val (sample_id), path("${sample_id}.g.vcf.gz")

    script:
    """
    samtools faidx ${reference}
    
    gatk CreateSequenceDictionary -R ${reference} -O ${reference.toString().replaceAll(/\.fa(sta)?$/, '.dict')}
    
    gatk --java-options "-Xmx4g" HaplotypeCaller \
    -R ${reference} \
    -I ${bam} \
    -O ${sample_id}.g.vcf.gz \
    -ERC GVCF
    """
}