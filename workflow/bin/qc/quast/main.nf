process QUAST {
    tag "QC_ASSEMBLE"
    
    publishDir "${params.outdir}/5-assemble/QUAST", mode: 'copy'
    
    input:
    tuple val(sample_id), path(contigs)
    tuple val(sample_id), path(scaffolds)
    tuple val(pair_id), path(trimmed_reads)

    output:
    tuple val(sample_id), path("quast_result_${sample_id}")

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
    ${scaffolds} 
    """
}