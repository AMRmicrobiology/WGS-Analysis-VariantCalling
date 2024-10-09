/*
DSL2 channels
*/
nextflow.enable.dsl=2

checkInputParams()

reference         = file("${params.reference}")

log.info """\

WGS - N F   P I P E L I N E - C L I N I C A L 
==============================================
Configuration environemnt:
    Out directory:             $params.outdir
    Fastq directory:           $params.input
    Reference directory:       $params.reference
"""
    .stripIndent()

//Call all the sub-work
include { FASTQC_QUALITY as FASTQC_QUALITY_ORIGINAL           }     from '../workflow/bin/qc/fastqc/main'
include { BUILD_INDEX as PERSONAL_GENOME_INDEX                }     from '../workflow/bin/bowtie/index/main'
include { BUILD_INDEX_2                                       }     from '../workflow/bin/bowtie/index/main_bwa'
include { TRIMMING                                            }     from '../workflow/bin/trimming/main'
include { FASTQC_QUALITY as FASTQC_QUALITY_FINAL              }     from '../workflow/bin/qc/fastqc/main'
include { ASSEMBLE                                            }     from '../workflow/bin/assemble/main'
include { PROKKA                                              }     from '../workflow/bin/anotations/prokka/main'
include { PERSONAL_GENOME_MAPPING                             }     from '../workflow/bin/bowtie/mapping/main'
include { BAKTA                                               }     from '../workflow/bin/anotations/bakta/main'
include { MARKDUPLICATE                                       }     from '../workflow/bin/gatk/picard/markduplicate/main'
include { ADDORREPLACE                                        }     from '../workflow/bin/gatk/picard/addorreplace/main'
include { HAPLOTYPECALLER                                     }     from '../workflow/bin/gatk/haplotype/main_1'
include { GENOTYPE as GENOTYPE_WILDTYPE                       }     from '../workflow/bin/gatk/genotype/main_1'
include { ALIGN as NORMALICE_WILDTYPE                         }     from '../workflow/bin/gatk/Filter/Align_1'
include { FILTER_VARIANTS as FILTER_VARIANTS_WILDTYPE         }     from '../workflow/bin/gatk/Filter/main_1'
include { ANOTATIONS as ANOTATION_SNPEFF                      }     from '../workflow/bin/snpeff/main'
include { AMR as POST-ANALYSIS-ABRICATE                       }     from '../workflow/bin/AMR/ABRIcate/main'
include { AMR_2 as POST-ANALYSIS-AMRFINDER                    }     from '../workflow/bin/AMR/AMRFinder/main'
include { AMR_3 as POST-ANALYSIS-DEEPARG                      }     from '../workflow/bin/AMR/DeepARG/main'


workflow {

//1st Step
//First Quality-control and build an INDEX - Specie Reference genome
    read_ch = Channel.fromFilePairs(params.input, size: 2 )
    fastqc_ch_original= FASTQC_QUALITY_ORIGINAL(read_ch.map{it -> it[1]})

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
    wildtype_only_ch = assemble_denovo_ch.contigs.first { it[0] ==~ /.*1$/ }

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
    gatk_haplotype_ch = HAPLOTYPECALLER (gatk_add_ch, params.personal_ref)

//GenotypeCaller 
//Perform joint genotyping 
    gatk_genotype_ch = GENOTYPE_WILDTYPE (gatk_haplotype_ch , params.personal_ref)


//Align
//This tool takes a VCF file, left-aligns the indels and trims common bases from indels, leaving them with a minimum representation.
//The same indel can often be placed at multiple positions and still represent the same haplotype.
//We are going to take the optionally splits multiallelic sites into biallelics and left-aligns individual alleles.
    aligns_and_normalized_ch = NORMALICE_WILDTYPE (gatk_genotype_ch, params.personal_ref)

//VatiantFilter
//Filter the VCF using the parametres to get a hight quality and cover in SNPs and INDELS "QUAL || MQ || DP ".
//all the parametres could be changen it, depends of the data.


    varaiant_filter_ch = FILTER_VARIANTS_WILDTYPE (aligns_and_normalized_ch, params.personal_ref)

//Genetic Variant annotation
//SNPEFF

    snpeff_anotation_ch = ANOTATION_SNPEFF (varaiant_filter_ch)

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