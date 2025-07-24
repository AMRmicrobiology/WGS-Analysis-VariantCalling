process ASSEMBLE {
    tag "Spades ${sample_id}"

    publishDir "${params.outdir}/assembly", mode: 'copy'
    input:

    tuple val (sample_id), path(pair_id_1), path(pair_id_2)

    output:

    tuple val (sample_id), path("${sample_id}.fasta"), emit: contigs
    tuple val (sample_id), path ("scaffolds_${sample_id}.fasta"), emit: scaffolds

    cpus 16
    memory '64 GB'
    time '24h'

    script:

    """
    spades.py -1 ${pair_id_1} -2 ${pair_id_2} --isolate -k auto -o ${sample_id}_spades_out && \
    mv ${sample_id}_spades_out/contigs.fasta ${sample_id}.fasta && \
    mv ${sample_id}_spades_out/scaffolds.fasta scaffolds_${sample_id}.fasta
    """
}
