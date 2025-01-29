process BUSCO {
    tag "GENOME COMPLETENESS ${sample_id}"

    container "$params.busco.docker"
    
    publishDir "${params.outdir}/BUSCO", mode: "copy" 

    input:
    tuple val(sample_id), path(assemble)

    output:
    tuple val(sample_id), path("${sample_code}_busco")

    script:

    """
    busco -i ${assemble} -m genome -o ${sample_id}_busco
    """
}