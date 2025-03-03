/*
DSL2 channels
*/
nextflow.enable.dsl=2

checkInputParams()

reference         = file("${params.reference}")


//Call all the sub-work
include { FASTQC_QUALITY as FASTQC_QUALITY_ORIGINAL           }     from '../bin/qc/fastqc/main'
include { TRIMMING                                            }     from '../bin/trimming/main'
include { FASTQC_QUALITY as FASTQC_QUALITY_FINAL              }     from '../bin/qc/fastqc/main'
include { MULTIQC                                             }     from '../bin/qc/multiqc/main'
include { BUILD_INDEX_1                                       }     from '../bin/bowtie/index/main_bwa'
include { BUILD_INDEX as PERSONAL_GENOME_INDEX                }     from '../bin/bowtie/index/main'
include { PERSONAL_GENOME_MAPPING                             }     from '../bin/bowtie/mapping/main'
include { MARKDUPLICATE                                       }     from '../bin/gatk/picard/markduplicate/main'
include { ADDORREPLACE                                        }     from '../bin/gatk/picard/addorreplace/main'
include { HAPLOTYPECALLER                                     }     from '../bin/gatk/haplotype/main_1'
include { GENOTYPE as GENOTYPE_ANALYSIS                       }     from '../bin/gatk/genotype/main'
include { ALIGN as NORMALICE_WILDTYPE                         }     from '../bin/gatk/Filter/align'
include { FILTER_VARIANTS as FILTER_VARIANTS_PARAM            }     from '../bin/gatk/Filter/main'
/*
include { ANOTATIONS as ANOTATION_SNPEFF                      }     from '../bin/snpeff/main'
*/

workflow reference {
    preprocess_output = workflow_pre_process()
    postprocess_output = workflow_post_process(preprocess_output.reference_ch, preprocess_output.fq_gz_reads_ch)
}

workflow workflow_pre_process {
    take:
    main:
    // Quality control and Index build
    read_ch = Channel.fromFilePairs(params.input, size: 2)
    /*
    fastqc_ch_original= FASTQC_QUALITY_ORIGINAL(read_ch.map{it -> it[1]})
    */
    // Trimming process
    trimmed_read_ch = TRIMMING(read_ch)
    fq_gz_reads_ch = trimmed_read_ch.trimmed_reads

    //Final Quality control after trimming
    fastq_ch_after = FASTQC_QUALITY_FINAL(trimmed_read_ch.trimmed_reads.map{it -> it[1]})
    
    //MULTIQC
    multiqc_ch = MULTIQC(fastqc_ch_original.qc_zip.collect(), fastq_ch_after.qc_zip.collect())

    //Reference Genome INDEX
    personal_ref_ch = Channel.fromPath(params.personal_ref)
    reference_ch = personal_ref_ch.map {
        ref -> 
        def sample_id = file(ref).baseName
        def ref_file = file(ref)
        return tuple (sample_id, ref_file)
    }

    personal_index_bwa_ch = BUILD_INDEX_1(reference_ch)
    personal_index_ch = PERSONAL_GENOME_INDEX(reference_ch)

    emit:
    reference_ch
    fq_gz_reads_ch
}

workflow workflow_post_process {

    take:
    reference_ch
    fq_gz_reads_ch
    
    main:

    //mapping process- Mapping used Specie ref. genome, include samtools sorted
    specie_mapping_ch = PERSONAL_GENOME_MAPPING(fq_gz_reads_ch, params.index_genome_personal)


    //Add groups and Mark duplicates
    bam_ch = specie_mapping_ch.map {
        tupla -> 
        def sample_id = tupla [0]
        def bam_path = tupla [2]
        return tuple (sample_id, bam_path)
    }

    gatk_mark_ch = MARKDUPLICATE (bam_ch)
    
    //Add or replace groups
    replace_ch = gatk_mark_ch.map {
        tupla -> 
        def sample_id = tupla [0]
        def replace_bam = tupla [1]
        return tuple (sample_id, replace_bam)
    }
    
    gatk_add_ch = ADDORREPLACE(replace_ch)

    //HAPLOTYPECALLER
    // realignment consistently incluide in the algoritme of GATK HaplotypeCaller.
    // minimum quality and confidence threshold are included
    haplotype_ch = gatk_add_ch
    .combine(reference_ch)
    .set { all_samples_ch }
    
    gatk_haplotype_ch = HAPLOTYPECALLER (all_samples_ch)

    //GenotypeCaller 
    //Perform joint genotyping 

    gatk_genotype_ch = GENOTYPE_ANALYSIS (gatk_haplotype_ch.out_files , gatk_haplotype_ch.reference_personal_genome)

    //Align
    //This tool takes a VCF file, left-aligns the indels and trims common bases from indels, leaving them with a minimum representation.
    //The same indel can often be placed at multiple positions and still represent the same haplotype.
    //We are going to take the optionally splits multiallelic sites into biallelics and left-aligns individual alleles.
    aligns_and_normalized_ch = NORMALICE_WILDTYPE (gatk_genotype_ch, gatk_haplotype_ch.reference_personal_genome)

    //VatiantFilter
    //Filter the VCF using the parametres to get a hight quality and cover in SNPs and INDELS "QUAL || MQ || DP ".
    //all the parametres could be changen it, depends of the data.
    varaiant_filter_ch = FILTER_VARIANTS_PARAM (aligns_and_normalized_ch, gatk_haplotype_ch.reference_personal_genome)
    /*
    // Decompress VCF
    vcf_ch = DECOMPRESS_VCF(variant_filter_ch.compl_vcf)

    
    // Functional annotation with SNPeff (optional)
    snpeff_ch = SNPEFF(agt_ch.combine_gff3, reference_ch, params.genome_name_db, agt_ch.protein_fasta, agt_ch.cds_fasta, vcf_ch)
    */
}

////////////////////////////////////////////////////////////////////////////////
// FUNCTIONS                                                                  //
////////////////////////////////////////////////////////////////////////////////


def checkInputParams() {
    // Check required parameters and display error messages
    boolean fatal_error = false
    if ( ! params.input) {
        log.warn("You need to provide a fastqDir (--fastqDir) or a bamDir (--bamDir)")
        fatal_error = true
    }
    if ( ! params.reference ) {
        log.warn("You need to provide a genome reference (--reference)")
        fatal_error = true
    }
    if (! params.personal_ref)  {
        log.warn("You need to provide a personal genome reference (--personal_ref)")
        fatal_error = true
    }
}