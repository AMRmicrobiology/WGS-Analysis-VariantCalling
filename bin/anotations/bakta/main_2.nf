process EXTRACT_CDS_FROM_BAKTA {
    tag "CDS EXTRACTOR for ${sample_id}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "docker://${params.agat.docker}" :
        params.agat.docker }"

    input:
    tuple val(sample_id), path(gff3_file), path(fna_file)

    output:
    path("cds_${sample_id}.fa"), emit: cds_fasta

    script:
    """
    sed 's/\\t?\\t/\\t.\\t/g' ${gff3_file} > cleaned_${sample_id}.gff3
    gffread cleaned_${sample_id}.gff3 -g ${fna_file} -x cds_${sample_id}.fa
    """
}
