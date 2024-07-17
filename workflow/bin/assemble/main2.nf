process ASSEMBLE_1 {
    tag "Sub-Sample-Spades ${sample_id}"

    publishDir "${params.outdir}/5-assemble", mode: 'copy'


    input:

    tuple val (sample_id), path(pair_id)

    output:

    tuple val (sample_id), path("${sample_id}_{1,2}.fastq")

    script:

    """
    seqtk sample -s100 ${pair_id[0]} 0.1 > ${sample_id}_1.fastq
    seqtk sample -s100 ${pair_id[1]} 0.1 > ${sample_id}_2.fastq
    """
}