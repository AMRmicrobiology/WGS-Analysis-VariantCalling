process FREEBAYES {
    tag "freebayes Variant Calling $sample_id"

    input:
    tuple val(sample_id), path(bam)
    path(reference)

    output:
    tuple val(sample_id), path("${sample_id}_freebayes.vcf")

    script:
    """
    # Ensure the reference is indexed
    samtools faidx ${reference}

    # Check if BAM is indexed, if not, index it
    if [ ! -f ${bam}.bai ]; then
        samtools index ${bam}
    fi

    # Run FreeBayes
    freebayes -f ${reference} ${bam} -v ${sample_id}_freebayes.vcf \\
    --ploidy 1 \\
    --min-coverage 10 \\
    --min-mapping-quality 20 \\
    --min-base-quality 20 \\
    --genotype-qualities \\
    --report-monomorphic \\
    --min-alternate-fraction 0.2 \\
    --use-mapping-quality \\
    --use-best-n-alleles 4
    """
}