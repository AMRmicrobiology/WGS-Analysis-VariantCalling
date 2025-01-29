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
    RGID= ${sample_id} \
    RGLB=lib1 \
    RGPL=ILLUMINA \
    RGPU=unit1 \
    RGSM=${sample_id} \
    CREATE_INDEX=True
    
    samtools flagstat ${sample_id}.RG.bam > ${sample_id}_samtools_flagstat.txt
    """
}