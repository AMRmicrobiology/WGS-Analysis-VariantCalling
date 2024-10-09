/*
DSL2 channels
*/
nextflow.enable.dsl=2

checkInputParams()

reference         = file("${params.reference}")

log.info """\

WGS - N F   P I P E L I N E
==============================================
Configuration environemnt:
    Out directory:             $params.outdir
    Fastq directory:           $params.input
    Reference directory:       $params.reference
"""
    .stripIndent()

//Call all the sub-work
include { FASTQC_QUALITY as FASTQC_QUALITY_ORIGINAL           }     from './workflow/bin/qc/fastqc/main'
include { BUILD_INDEX as PERSONAL_GENOME_INDEX                }     from './workflow/bin/bowtie/index/main'
include { BUILD_INDEX_2                                       }     from './workflow/bin/bowtie/index/main_bwa'
include { TRIMMING                                            }     from './workflow/bin/trimming/main'
include { FASTQC_QUALITY as FASTQC_QUALITY_FINAL              }     from './workflow/bin/qc/fastqc/main'
include { ASSEMBLE                                            }     from './workflow/bin/assemble/main'
include { PERSONAL_GENOME_MAPPING                             }     from './workflow/bin/bowtie/mapping/main'
include { MARKDUPLICATE                                       }     from './workflow/bin/gatk/picard/markduplicate/main'
include { ADDORREPLACE                                        }     from './workflow/bin/gatk/picard/addorreplace/main'
include { HAPLOTYPECALLER                                     }     from './workflow/bin/gatk/haplotype/main'
include { COMBINE_GVCFS as COMBINE_GVCFS_CP14                 }     from './workflow/bin/gatk/combine/main'
include { COMBINE_GVCFS as COMBINE_GVCFS_NET                  }     from './workflow/bin/gatk/combine/main'
include { GENOTYPE as GENOTYPE_CP14                           }     from './workflow/bin/gatk/genotype/main'
include { GENOTYPE as GENOTYPE_NET                            }     from './workflow/bin/gatk/genotype/main'
include { GENOTYPE as GENOTYPE_WILDTYPE                       }     from './workflow/bin/gatk/genotype/main'
include { ALIGN as NORMALICE_CP14                             }     from './workflow/bin/gatk/Filter/Align'
include { ALIGN as NORMALICE_NET                              }     from './workflow/bin/gatk/Filter/Align'
include { ALIGN as NORMALICE_WILDTYPE                         }     from './workflow/bin/gatk/Filter/Align'
include { FILTER_VARIANTS as FILTER_VARIANTS_CP14             }     from './workflow/bin/gatk/Filter/main'
include { FILTER_VARIANTS as FILTER_VARIANTS_NET              }     from './workflow/bin/gatk/Filter/main'
include { FILTER_VARIANTS as FILTER_VARIANTS_WILDTYPE         }     from './workflow/bin/gatk/Filter/main'
include { ANOTATIONS as ANOTATION_SNPEFF                      }     from './workflow/bin/snpeff/main'
include { AMR as POST-ANALYSIS-ABRICATE                       }     from './workflow/bin/AMR/ABRIcate/main'
include { AMR_2 as POST-ANALYSIS-AMRFINDER                    }     from './workflow/bin/AMR/AMRFinder/main'
include { AMR_3 as POST-ANALYSIS-DEEPARG                      }     from './workflow/bin/AMR/DeepARG/main'


workflow {

//1st Step
//First Quality-control and build an INDEX - Specie Reference genome
    read_ch = Channel.fromFilePairs(params.input, size: 2 )
/*
    fastqc_ch_original= FASTQC_QUALITY_ORIGINAL(read_ch.map{it -> it[1]})
*/

//Pruning (Bowtie2 index + Trimming)
    //Trimming-Reads - Cleaning paired reads and trimming adapters
    trimmed_read_ch = TRIMMING(read_ch)
    //Build an INDEX - personal reference genome
    personal_ref_ch = Channel.fromPath( [ "$params.personal_ref" ] )
    personal_index_ch  = PERSONAL_GENOME_INDEX (personal_ref_ch)
    personal_index_bwa_ch = BUILD_INDEX_2 (params.personal_ref)

    //Final Quality control after trimming
    FASTQC_QUALITY_FINAL(trimmed_read_ch.trimmed_reads.map{it -> it[1]})

// DE NOVO ASSEMBLE
    assemble_denovo_ch = ASSEMBLE(trimmed_read_ch.trimmed_reads)

//2nd Step
//mapping process- Mapping used Specie ref. genome, include samtools sorted
    specie_mapping_ch   = PERSONAL_GENOME_MAPPING(trimmed_read_ch.trimmed_reads, params.index_genome_personal)
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
    gatk_wild_type_ch = HAPLOTYPECALLER (gatk_add_ch, params.personal_ref)


    gvcf_CP14_ch = gatk_wild_type_ch.filter { it.name =~ /AB([2-9]|10|11)\.g\.vcf\.gz/ }.collect()
    gvcf_NET_ch = gatk_wild_type_ch.filter { it.name =~ /AB(1[2-9]|20|21)\.g\.vcf\.gz/ }.collect()
    gvcf_WILTYPE_ch = gatk_wild_type_ch.filter {it.name =~ /AB([1])\.g\.vcf\.gz/}

//CombineGVCF
//Combine per-sample gVCF files produced by HaplotypeCaller into a multi-sample gVCF file.
//The main advantage is the ability to combine multiple intervals at once without building a GenomicsDB. 

    combine_C14_ch= COMBINE_GVCFS_CP14 (gvcf_CP14_ch, params.personal_ref, 'CP14')
    combine_NET_ch= COMBINE_GVCFS_NET (gvcf_NET_ch, params.personal_ref, 'NET')

//GenotypeCaller 
//Perform joint genotyping 
    genotype_CP14_ch = GENOTYPE_CP14 (combine_C14_ch , params.personal_ref, 'CP14')
    genotype_NET_ch = GENOTYPE_NET (combine_NET_ch, params.personal_ref, 'NET')
    gentoype_wildtype_ch = GENOTYPE_WILDTYPE (gvcf_WILTYPE_ch, params.personal_ref, 'WILTYPE')


//Align
//This tool takes a VCF file, left-aligns the indels and trims common bases from indels, leaving them with a minimum representation.
//The same indel can often be placed at multiple positions and still represent the same haplotype.
//We are going to take the optionally splits multiallelic sites into biallelics and left-aligns individual alleles.
    aligns_and_normalized_CP14_ch = NORMALICE_CP14 (genotype_CP14_ch, params.personal_ref, 'CP14')
    aligns_and_normalized_NET_ch = NORMALICE_NET (genotype_NET_ch, params.personal_ref, 'NET')
    aligns_and_normalized_WILDTYPE_ch = NORMALICE_WILDTYPE (gentoype_wildtype_ch, params.personal_ref, 'WILDTYPE')

//VatiantFilter
//Filter the VCF using the parametres to get a hight quality and cover in SNPs and INDELS "QUAL || MQ || DP ".
//all the parametres could be changen it, depends of the data.

    varaiant_filter_CP14_ch = FILTER_VARIANTS_CP14 (aligns_and_normalized_CP14_ch, params.personal_ref, 'CP14')
    varaiant_filter_NET_ch = FILTER_VARIANTS_NET (aligns_and_normalized_NET_ch, params.personal_ref, 'NET')
    varaiant_filter_WILDTYPE_ch = FILTER_VARIANTS_WILDTYPE (aligns_and_normalized_WILDTYPE_ch, params.personal_ref, 'WILDTYPE')

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