process QUAST {
    tag "QC_ASSEMBLE"
    
    input:
    tuple val(sample_id), path(contigs), path(scaffolds)
    tuple val(pair_id), path(trimmed_reads)

    output:
    tuple val(sample_id), path("quast_result_${sample_id}")

    script:
    """
    quast.py \\
    -o quast_result_${sample_id} \\
    -m 0 -t 2 -k \\
    --k-mer-size 127 \\
    --circos \\
    --pe1 ${trimmed_reads[0]} \\
    --pe2 ${trimmed_reads[1]} \\
    ${contigs} \\
    ${scaffolds} 
    """
}
