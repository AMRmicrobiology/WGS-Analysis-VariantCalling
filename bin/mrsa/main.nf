process MRSA {
    tag "MRSA process SPATYPER-SCCMEC ${sample_id}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "docker://${params.short_wgs.docker}" :
        params.short_wgs.docker }"

    publishDir "${params.outdir}/5-MRSA/spaTyper" , mode:"copy"

    input:

    tuple val (sample_id), path(contigs)


    output:

    tuple val (sample_id), path ("${sample_id}_spatype.txt")


    script:
    
    """
    download-spatypes.sh
    
    spaTyper -d /opt/conda/envs/env/share/spatyper-0.3.3 -f ${contigs} --output ${sample_id}_spatype.txt 

    """
}

process SCCMEC {
    tag "MRSA process SPATYPER-SCCMEC ${sample_id}"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "docker://${params.short_wgs.docker}" :
        params.short_wgs.docker }"

    publishDir "${params.outdir}/5-MRSA/SCCMEC" , mode: "copy"
    
    input:

    tuple val (sample_id), path(contigs)


    output:

    path "*.tsv"


    script:
    
    """
    sccmec --input ${contigs} --prefix ${sample_id}

    """
}