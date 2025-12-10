process ASSEMBLE {
    tag "Spades assembly for ${sample_id}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "docker://${params.short_wgs.docker}" :
        params.short_wgs.docker }"
    
    input:

    tuple val (sample_id), path(pair_id)

    output:

    tuple val (sample_id), path("${sample_id}.fasta"), emit: contigs
    tuple val (sample_id), path ("scaffolds_${sample_id}.fasta"), emit: scaffolds

    cpus 16
    memory '64 GB'
    time '24h'

    script:

    """
    spades.py -1 ${pair_id[0]} -2 ${pair_id[1]} --isolate -k auto -o ${sample_id}_spades_out && \
    mv ${sample_id}_spades_out/contigs.fasta ${sample_id}.fasta && \
    mv ${sample_id}_spades_out/scaffolds.fasta scaffolds_${sample_id}.fasta
    """
}
