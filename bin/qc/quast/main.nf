process QUAST {
    tag "QC_ASSEMBLE"
    
    publishDir "${params.outdir}/1-QC/genomeQC/QUAST", mode: 'copy'
    
    input:
    tuple val(sample_id), path(contigs), path(scaffolds), path(trimmed_reads)

    output:
    tuple val(sample_id), path("quast_result_${sample_id}/report_${sample_id}.tsv"), emit: report_tsv_quast
    tuple val(sample_id), path("quast_result_${sample_id}/report.txt"), emit: report_txt_quast
    tuple val(sample_id), path("quast_result_${sample_id}"), emit: direct_quast    

    script:

    """
    quast.py \\
    -o quast_result_${sample_id} \\
    -m 500 -t 4 -k \\
    --k-mer-size 127 \\
    --circos \\
    --pe1 ${trimmed_reads[0]} \\
    --pe2 ${trimmed_reads[1]} \\
    --gene-finding \\
    --rna-finding \\
    --contig-thresholds 0 \\
    ${contigs} \\
    ${scaffolds} \\
    --threads 8

    mv quast_result_${sample_id}/report.tsv quast_result_${sample_id}/report_${sample_id}.tsv
    """
}