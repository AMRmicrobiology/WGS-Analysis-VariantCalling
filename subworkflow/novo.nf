/*
DSL2 channels
*/
nextflow.enable.dsl=2

checkInputParams()

reference         = file("${params.reference}")

log.info """\

WGS - P A R A M E T R E S
==============================================
Configuration environemnt:
    Out directory:             $params.outdir
    Fastq directory:           $params.input
    Reference directory:       $params.reference
"""
    .stripIndent()

//Call all the sub-work

include { FASTQC_QUALITY as FASTQC_QUALITY_ORIGINAL           }     from '../workflow/bin/qc/fastqc/main'
include { TRIMMING                                            }     from '../workflow/bin/trimming/main'
include { FASTQC_QUALITY as FASTQC_QUALITY_FINAL              }     from '../workflow/bin/qc/fastqc/main'
include { ASSEMBLE                                            }     from '../workflow/bin/assemble/main'
include { QUAST                                               }     from '../workflow/bin/qc/quast/main'
include { MULTIQC                                             }     from '../workflow/bin/qc/multiqc/main' 
include { PROKKA                                              }     from '../workflow/bin/anotations/prokka/main'
include { BAKTA                                               }     from '../workflow/bin/anotations/bakta/main'
include { BUILD_INDEX_1                                       }     from '../workflow/bin/bowtie/index/main_bwa'
include { BUILD_INDEX as PERSONAL_GENOME_INDEX                }     from '../workflow/bin/bowtie/index/main'
include { PERSONAL_GENOME_MAPPING                             }     from '../workflow/bin/bowtie/mapping/main'
include { MARKDUPLICATE                                       }     from '../workflow/bin/gatk/picard/markduplicate/main'
include { ADDORREPLACE                                        }     from '../workflow/bin/gatk/picard/addorreplace/main'
include { HAPLOTYPECALLER                                     }     from '../workflow/bin/gatk/haplotype/main_1'
include { GENOTYPE as GENOTYPE_ANALYSIS                       }     from '../workflow/bin/gatk/genotype/main_1'
include { ALIGN as NORMALICE_WILDTYPE                         }     from '../workflow/bin/gatk/Filter/Align_1'
include { FILTER_VARIANTS as FILTER_VARIANTS_PARAM            }     from '../workflow/bin/gatk/Filter/main_1'
include { AGT                                                 }     from '../workflow/bin/anotations/main'
include { DECOMPRESS_VCF                                      }     from '../workflow/bin/snpeff/main_2'
include { SNPEFF                                              }     from '../workflow/bin/snpeff/main'
include { AMR as POST_ANALYSIS_ABRICATE                       }     from '../workflow/bin/AMR/abricate/main'
include { AMR_2 as POST_ANALYSIS_AMRFINDER                    }     from '../workflow/bin/AMR/AMRFinder/main'
include { AMR_3 as POST_ANALYSIS_DEEPARG                      }     from '../workflow/bin/AMR/DeepARG/main'
/*



*/
workflow novo {
    preprocess_output = workflow_pre_process()
    postprocess_output = workflow_post_process(preprocess_output.personal_ref_ch, preprocess_output.fq_gz_reads_ch)
    amrprocess_output = workflow_amr( preprocess_output.contigs_ch)
}

workflow workflow_pre_process {

    take:
    main:
    // Quality control y construcción del índice
    read_ch = Channel.fromFilePairs(params.input, size: 2)

    fastqc_ch_original= FASTQC_QUALITY_ORIGINAL(read_ch.map{it -> it[1]})

    // Trimming de las lecturas
    trimmed_read_ch = TRIMMING(read_ch)
    fq_gz_reads_ch = trimmed_read_ch.trimmed_reads
   
    //Final Quality control after trimming
    fastq_ch_after = FASTQC_QUALITY_FINAL(trimmed_read_ch.trimmed_reads.map{it -> it[1]})

    // Ensamblado de novo
    assemble_denovo_ch = ASSEMBLE(trimmed_read_ch.trimmed_reads)
    wildtype_only_ch = assemble_denovo_ch.contigs.first { it[0] ==~ /.*[^0-9]1$/ }
    contigs_ch = assemble_denovo_ch.contigs
    scaffolds_ch = assemble_denovo_ch.scaffolds

    //QUAST
    quast_ch = QUAST(assemble_denovo_ch.contigs, assemble_denovo_ch.scaffolds, trimmed_read_ch.trimmed_reads )
    .map { tuple -> tuple[1] } 
    //MULTIQC
    multiqc_ch = MULTIQC(fastqc_ch_original.qc_zip.collect(), fastq_ch_after.qc_zip.collect(), quast_ch.collect())

    // Construir un índice del genoma de referencia personal
    personal_ref_ch = wildtype_only_ch
    personal_index_bwa_ch = BUILD_INDEX_1(personal_ref_ch)
    personal_index_ch = PERSONAL_GENOME_INDEX(personal_ref_ch)

    emit:
    contigs_ch
    scaffolds_ch
    personal_ref_ch
    fq_gz_reads_ch
}

workflow workflow_post_process {
    
    take:
    personal_ref_ch
    fq_gz_reads_ch
    
    main:
    //anotations process
    prokka_anotation_ch = PROKKA(personal_ref_ch)
    //anotations BAKTA
    bakta_anotation_ch = BAKTA(personal_ref_ch)
    //merge anotations
    agt_ch = AGT(prokka_anotation_ch.prokka_gff, bakta_anotation_ch.bakta_gff3, personal_ref_ch)

    //2nd Step
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
    gatk_haplotype_ch = HAPLOTYPECALLER (gatk_add_ch, personal_ref_ch)

    //GenotypeCaller 
    //Perform joint genotyping 
    gatk_genotype_ch = GENOTYPE_ANALYSIS (gatk_haplotype_ch , personal_ref_ch)

    //Align
    //This tool takes a VCF file, left-aligns the indels and trims common bases from indels, leaving them with a minimum representation.
    //The same indel can often be placed at multiple positions and still represent the same haplotype.
    //We are going to take the optionally splits multiallelic sites into biallelics and left-aligns individual alleles.
    aligns_and_normalized_ch = NORMALICE_WILDTYPE (gatk_genotype_ch, personal_ref_ch)

    //VatiantFilter
    //Filter the VCF using the parametres to get a hight quality and cover in SNPs and INDELS "QUAL || MQ || DP ".
    //all the parametres could be changen it, depends of the data.
    varaiant_filter_ch = FILTER_VARIANTS_PARAM (aligns_and_normalized_ch, personal_ref_ch)

    //DESCROMPRES VCF
    vcf_ch = DECOMPRESS_VCF(varaiant_filter_ch.compl_vcf)

    //SNPeFF
    //Funcional anotations
    snpeff_ch = SNPEFF(agt_ch.combine_gff3, personal_ref_ch, params.genome_name_db, agt_ch.protein_fasta, agt_ch.cds_fasta, vcf_ch)
    
}

workflow workflow_amr {
    take:
    contigs_ch
    
    main:
    //AMR
    //AMR1-ABRIcate
    abricate_ch = POST_ANALYSIS_ABRICATE(contigs_ch)

    //AMR2-RESFINDER
    resfinder_ch = POST_ANALYSIS_AMRFINDER(contigs_ch)

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