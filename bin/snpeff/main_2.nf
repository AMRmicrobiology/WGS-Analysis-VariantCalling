process DECOMPRESS_VCF {
    tag "DESCOMPRESS VCF FOR SNPEFF"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "docker://${params.short_wgs.docker}" :
        params.short_wgs.docker }"
        
    input:
    tuple val (sample_id), path(filtered_vcf)

    output:
    tuple val(sample_id), path("${filtered_vcf.baseName}")
    
    script:
    """
    gunzip -c ${filtered_vcf} > ${filtered_vcf.baseName}
    """
}