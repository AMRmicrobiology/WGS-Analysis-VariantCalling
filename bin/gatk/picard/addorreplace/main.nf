process ADDORREPLACE {
    tag "Add or Replace $sample_id "

    input:
    tuple val(sample_id), path(replace)

    output:
    tuple val(sample_id), path("${sample_id}.RG.bam")
    

    script:

    """
    picard AddOrReplaceReadGroups \
    INPUT=${replace} \
    OUTPUT=${sample_id}.RG.bam \
    RGID= M06646 \
    RGLB=lib1 \
    RGPL=ILLUMINA \
    RGPU=unit1 \
    RGSM=${sample_id} \
    CREATE_INDEX=True
    """
}