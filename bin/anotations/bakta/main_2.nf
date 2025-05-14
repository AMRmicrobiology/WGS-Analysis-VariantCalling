process EXTRACT_CDS_FROM_BAKTA {
    tag "CDS EXTRACTOR for ${sample_id}"

    container "$params.agat.docker"


    input:
    tuple val(sample_id), path(gff3_file), path(fna_file)

    output:
    tuple val(sample_id), path("cds_${sample_id}.fa"), emit: cds_fasta

    script:
    """
    gffread ${gff3_file} -g ${fna_file} -x cds_${sample_id}.fa
    """

}
