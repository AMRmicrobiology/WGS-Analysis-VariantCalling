process TRIMMING {
    tag "fastp ${pair_id}"

    input:
    tuple val(pair_id), path (reads) 

    output:
    tuple val(pair_id), path("${pair_id}_clean_{1,2}.fq.gz"), emit: trimmed_reads
    path("*.html"), emit: report

    script: 
    """
    fastp\\
        -i ${reads[0]}\\
        -I ${reads[1]}\\
        -o ${pair_id}_clean_1.fq.gz\\
        -O ${pair_id}_clean_2.fq.gz\\
        --trim_poly_g\\
        --trim_poly_x\\
        --detect_adapter_for_pe\\
        --cut_front 15\\
        --cut_tail 20\\
        --cut_mean_quality 20\\
        --length_required 50\\
        -h out_${pair_id}_Fastp.html
    """
}