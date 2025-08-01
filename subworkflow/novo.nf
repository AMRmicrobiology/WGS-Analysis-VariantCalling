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
    DB SNPeFF name:            $params.genome_name_db
"""
    .stripIndent()

//Call all the sub-work

include { FASTQC_QUALITY as FASTQC_QUALITY_ORIGINAL           }     from '../bin/qc/fastqc/main'
include { TRIMMING                                            }     from '../bin/trimming/main'
include { FASTQC_QUALITY as FASTQC_QUALITY_FINAL              }     from '../bin/qc/fastqc/main'
include { MULTIQC                                             }     from '../bin/qc/multiqc/main' 
include { KRAKEN;SEQTK_PRUNE                                  }     from '../bin/kraken/main'
include { ASSEMBLE                                            }     from '../bin/assemble/main'
include { FILTER_CONTIGS                                      }     from '../bin/qc/polish/filter'
include { ALIGMENT_PILON;PILON_POLISH                         }     from '../bin/qc/polish/main'
include { PROKKA                                              }     from '../bin/anotations/prokka/main'
include { BAKTA                                               }     from '../bin/anotations/bakta/main'
include { QUAST                                               }     from '../bin/qc/quast/main'
include { BUSCO                                               }     from '../bin/qc/busco/main'
include { MULTIQC_2 as POST_MULTIQC                           }     from '../bin/qc/multiqc/main_2' 
include { BUILD_INDEX_1                                       }     from '../bin/bowtie/index/main_bwa'
include { BUILD_INDEX as PERSONAL_GENOME_INDEX                }     from '../bin/bowtie/index/main'
include { AGT                                                 }     from '../bin/anotations/main'
include { PERSONAL_GENOME_MAPPING                             }     from '../bin/bowtie/mapping/main'
include { MARKDUPLICATE                                       }     from '../bin/gatk/picard/markduplicate/main'
include { ADDORREPLACE                                        }     from '../bin/gatk/picard/addorreplace/main'
include { HAPLOTYPECALLER                                     }     from '../bin/gatk/haplotype/main_3'
include { GENOTYPE as GENOTYPE_ANALYSIS                       }     from '../bin/gatk/genotype/main'
include { ALIGN as NORMALICE_WILDTYPE                         }     from '../bin/gatk/Filter/align'
include { FILTER_VARIANTS as FILTER_VARIANTS_PARAM            }     from '../bin/gatk/Filter/main'
include { DECOMPRESS_VCF                                      }     from '../bin/snpeff/main_2'
include { SNPEFF                                              }     from '../bin/snpeff/main'
include { AMR as POST_ANALYSIS_ABRICATE                       }     from '../bin/AMR/abricate/main'
include { AMR_2 as POST_ANALYSIS_AMRFINDER                    }     from '../bin/AMR/AMRFinder/main'


workflow novo {
    preprocess_output = workflow_pre_process()
    anotationprocess_output = workflow_anotation_process(preprocess_output.personal_ref_ch, preprocess_output.accurance_fasta_ch, preprocess_output.wildtype_only_ch)
    mappingprocess_output = workflow_mapping_process(preprocess_output.fq_gz_reads_ch, preprocess_output.personal_ref_ch, anotationprocess_output.accurance_fasta_ch)
    /*
    amrprocess_output = workflow_amr( preprocess_output.contigs_ch)*/
}

workflow workflow_pre_process {

    main:
    // Quality control and index build
    read_ch = Channel.fromFilePairs(params.input, size: 2)

    fastqc_ch_original= FASTQC_QUALITY_ORIGINAL(read_ch.map{it -> it[1]})

    // Trimming process
    trimmed_read_ch = TRIMMING(read_ch)
    fq_gz_reads_ch = trimmed_read_ch.trimmed_reads

    //KRAKEN
    kraken_ch = KRAKEN(fq_gz_reads_ch)

    //Final Quality control after trimming
    fastq_ch_after = FASTQC_QUALITY_FINAL(trimmed_read_ch.trimmed_reads.map{it -> it[1]})

    //PRUNNING
    fastq_prunning_ch = fq_gz_reads_ch.join(kraken_ch.keep_ids).map {
        sample_id,reads_pair, keep_ids ->
        def (r1, r2) = reads_pair
        tuple (sample_id, [r1, r2], keep_ids)
    }
    
    prune_ch = SEQTK_PRUNE(fastq_prunning_ch)
   
    //de novo assemble
    assemble_denovo_ch = ASSEMBLE(prune_ch)
    contigs_ch = assemble_denovo_ch.contigs
    scaffolds_ch = assemble_denovo_ch.scaffolds
    
    //Filter seq low quality contigs
    filtered_contigs_ch = FILTER_CONTIGS(contigs_ch)
 
    //Polishing Illumina SEQ
    polish_data_ch = filtered_contigs_ch
        .join(trimmed_read_ch.trimmed_reads)
        .map { sample_id, contigs, reads_clean_pair -> 
        def (r1, r2) = reads_clean_pair
        tuple (sample_id, contigs , [r1, r2])
    }

    polishing_illumina_ch = ALIGMENT_PILON(polish_data_ch)
    
    polish_data_index_ch = filtered_contigs_ch
        .join(polishing_illumina_ch.aln_bam)
        .map { sample_id, contigs, index_bam -> 
        tuple (sample_id, contigs , index_bam)
    }

    pilon_polish_ch = PILON_POLISH(polish_data_index_ch)
    accurance_fasta_ch = pilon_polish_ch.pilon_fa

    wildtype_only_ch = accurance_fasta_ch.first { it[0] == params.wildtype_code }

    // Index build
    personal_ref_ch = wildtype_only_ch
    personal_index_bwa_ch = BUILD_INDEX_1(personal_ref_ch)
    personal_index_ch = PERSONAL_GENOME_INDEX(personal_ref_ch)

    //Emit results
    emit:
    scaffolds_ch
    accurance_fasta_ch
    fq_gz_reads_ch
    personal_ref_ch
    wildtype_only_ch

}


workflow workflow_anotation_process {

    take:
    accurance_fasta_ch
    personal_ref_ch
    wildtype_only_ch


    main:

    //PROKKA
    prokka_annotation_ch = PROKKA(wildtype_only_ch)
    bakta_annotation_ch = BAKTA(wildtype_only_ch)
    
    //merge anotations
    agt_ch = AGT(prokka_annotation_ch.prokka_gff, bakta_annotation_ch.bakta_gff3, wildtype_only_ch)

}

workflow workflow_mapping_process {

    take:
    fq_gz_reads_ch
    personal_ref_ch
    accurance_fasta_ch

    main:

    //mapping process- Mapping used Specie ref. genome, include samtools sorted
    specie_mapping_ch = PERSONAL_GENOME_MAPPING(fq_gz_reads_ch, params.index_genome_personal)

    //Add groups and Mark duplicates

    bam_ch = specie_mapping_ch.all_outputs.map { sample_id, sam, bam, bai, metrics, flagstat -> 
    tuple(sample_id, bam)
    }
    //Add or replace read groups
    gatk_add_ch = ADDORREPLACE(bam_ch)
    
    //Mark duplicates
    gatk_mark_ch = MARKDUPLICATE (gatk_add_ch)
   
    //GATK PROCESS
    //GATK PREPARE CHANNEL
    replace_ch = gatk_mark_ch.dedup_bam
    //HAPLOTYPECALLER
    // realignment consistently incluide in the algoritme of GATK HaplotypeCaller.
    // minimum quality and confidence threshold are included
    gatk_haplotype_ch = HAPLOTYPECALLER (replace_ch, personal_ref_ch)
   
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