process QUAST {
    tag "QC_ASSEMBLE"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "docker://${params.short_wgs.docker}" :
        params.short_wgs.docker }"
    
    publishDir "${params.outdir}/1-QC/genomeQC/QUAST", mode: 'copy'
    
    errorStrategy 'ignore'
    
    input:
    tuple val(sample_id), path(contigs), path(trimmed_reads)

    output:

    tuple val(sample_id), path("quast_result_${sample_id}/") 

    script:

    """
    quast.py \\
    -o quast_result_${sample_id} \\
    -m 500 \\
    --threads 8 \\
    --k-mer-size 127 \\
    --circos \\
    --pe1 ${trimmed_reads[0]} \\
    --pe2 ${trimmed_reads[1]} \\
    --gene-finding \\
    --rna-finding \\
    --contig-thresholds 0 \\
    ${contigs} \\

    """
}