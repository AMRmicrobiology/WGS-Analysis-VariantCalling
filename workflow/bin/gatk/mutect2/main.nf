process MUTECT2 {
    tag ""


    input:



    output:



    script:
    """
    Mutect2 \
    -R reference.fa \
    -I sample.bam \
    -O sample.vcf.gz \
    --max-reads-per-alignment-start 0 \
    --min-base-quality-score 20 \
    --initial-tumor-lod 0 \
    --tumor-lod-to-emit 0 \
    --af-of-alleles-not-in-resource 4e-3 \
    --pruning-lod-threshold -4
    
    """
}