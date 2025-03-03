process MRSA {
    tag "MRSA process SPATYPER-SCCMEC ${sample_id}"

    


    input:

    tuple val (sample_id), path(contings)


    output:

    tuple val (sample_id), path ()


    script:
    
    """
    
    sccmec --input ${contings} --prefix ${sample_id}
    
    
    """
}