process HAPLOTYPECALLER {
    tag "Haplotype ${sample_id}"
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "docker://${params.gatk4.docker}" :
        params.gatk4.docker }"

    publishDir "${params.outdir}", mode: 'copy', saveAs: { filename ->
        if (filename.endsWith(".vcf.")) "8-tryVCF/VCF/$filename"
        else null
    }

    input:
    tuple val (sample_id), path (bam), val(wt), path (reference)

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