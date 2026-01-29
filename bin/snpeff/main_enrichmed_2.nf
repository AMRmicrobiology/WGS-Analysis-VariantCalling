process JOIN_SNPEFF_GFF {
  tag "${sample_id}"

  publishDir "${params.outdir}/Variant_annotations_enrichment", mode: 'copy'

  input:
    tuple val(sample_id), path(vcf), path(gff)

  output:
    tuple val(sample_id), path("${sample_id}.snpeff.prokka.tsv"), path("${sample_id}.snpeff.prokka.report.txt")

  script:
    """
    python3 ${workflow.projectDir}/bin/snpeff/vcf_snpeff_gff_join.py \
      --vcf ${vcf} \
      --gff ${gff} \
      --out ${sample_id}.snpeff.prokka.tsv
    """
}