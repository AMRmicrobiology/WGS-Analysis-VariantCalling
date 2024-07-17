process ASSEMBLE {
    tag "Spades ${sample_id}"

    publishDir "${params.outdir}/5-assemble", mode: 'copy'


    input:

    tuple val (sample_id), path(pair_id)

    output:

    tuple val (sample_id), path("${sample_id}.fa"), path ("scaffolds_${sample_id}.fasta"), emit: contings_scaffolds

    cpus 16
    memory '64 GB'
    time '24h'

    script:

    """
    spades.py -1 ${pair_id[0]} -2 ${pair_id[1]} --isolate -k auto -o ${sample_id}_spades_out
    mv ${sample_id}_spades_out/contigs.fasta ${sample_id}.fa
    mv ${sample_id}_spades_out/scaffolds.fasta scaffolds_${sample_id}.fasta
    """
}