process DECOMPRESS_VCF {
    tag "DESCOMPRESS VCF FOR SNPEFF"
    input:
    tuple val (sample_id), path(filtered_vcf)

    output:
    tuple val(sample_id), path("${filtered_vcf.baseName}")
    script:
    """
    # Descomprimir el archivo VCF.gz
    gunzip -c ${filtered_vcf} > ${filtered_vcf.baseName}
    """
}