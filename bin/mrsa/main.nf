process MRSA {
    tag "MRSA process SPATYPER-SCCMEC ${sample_id}"


    input:

    tuple val (sample_id), path(contigs)


    output:

    tuple val (sample_id), path ()


    script:
    
    """

    // SE LE PUEDE ANADIR UN COLLECT

    download-spatypes.sh
    
    spaTyper -d /opt/conda/envs/env/share/spatyper-0.3.3 -f ${contigs} --output spatype.txt 

    """
}

process SCCMEC {
    tag "MRSA process SPATYPER-SCCMEC ${sample_id}"

    
    input:

    tuple val (sample_id), path(contigs)


    output:

    tuple val (sample_id), path ()


    script:
    
    """   
    sccmec --input ${contings} --prefix ${sample_id}

    """
}