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
    QCprocess_output = workflow_QC_process(preprocess_output.busco_ch, preprocess_output.quast_ch)
    postprocess_output = workflow_post_process(preprocess_output.personal_ref_ch, preprocess_output.fq_gz_reads_ch, preprocess_output.accurance_fasta_ch)
    /*amrprocess_output = workflow_amr( preprocess_output.contigs_ch)*/
}

workflow workflow_pre_process {

    take:
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

    //PROKKA
    prokka_annotation_ch = PROKKA(accurance_fasta_ch)
    prokka_fna_ch = accurance_fasta_ch.join(prokka_annotation_ch.prokka_path).map {
        sample_id, contigs, prokka_fna -> 
        tuple (sample_id, prokka_fna)
    }

    wildtype_only_ch = prokka_annotation_ch.prokka_path.filter { it[0] == params.wildtype_code }

    //BAKTA
    bakta_annotation_ch = BAKTA(prokka_fna_ch)

    //BUSCO
    busco_ch = BUSCO(accurance_fasta_ch)

    //QUAST

    quast_input_ch = accurance_fasta_ch.join(trimmed_read_ch.trimmed_reads)
        .map { sample_id, contigs, reads_clean_pair ->
        def (r1, r2) = reads_clean_pair
        tuple (sample_id, contigs, [r1, r2])
    }

    quast_ch = QUAST(quast_input_ch)

    //MULTIQC
    multiqc_ch = MULTIQC(fastqc_ch_original.qc_zip.collect(), fastq_ch_after.qc_zip.collect())
    
    // Index build
    personal_ref_ch = wildtype_only_ch
    personal_index_bwa_ch = BUILD_INDEX_1(personal_ref_ch)
    personal_index_ch = PERSONAL_GENOME_INDEX(personal_ref_ch)

    //merge anotations
    wildtype_prokka_ch = prokka_annotation_ch.prokka_gff.filter { it.endsWith("${params.wildtype_code}.gff") }
    wildtype_bakta_ch = bakta_annotation_ch.bakta_gff3.filter { it.endsWith("${params.wildtype_code}.gff3") }
    
    agt_ch = AGT(wildtype_prokka_ch, wildtype_bakta_ch, personal_ref_ch)

    //Emit results
    emit:
    contigs_ch
    scaffolds_ch
    accurance_fasta_ch
    fq_gz_reads_ch
    busco_ch
    quast_ch
    personal_ref_ch

}

workflow workflow_QC_process {

    take:
    busco_ch
    quast_ch

    main:
    multiqc_2_ch = POST_MULTIQC(quast_ch.map{ it -> it[1] }.collect(), busco_ch.map{ it -> it[1] }.collect())

}


workflow workflow_post_process {
    
    take:
    accurance_fasta_ch
    fq_gz_reads_ch
    personal_ref_ch
        
    main:
    //Mapping process- Mapping used Specie ref. genome, include samtools sorted
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
/*
    //DESCROMPRES VCF
    vcf_ch = DECOMPRESS_VCF(varaiant_filter_ch.compl_vcf)

    //SNPeFF
    //Funcional anotations
    snpeff_ch = SNPEFF(agt_ch.combine_gff3, personal_ref_ch, params.genome_name_db, agt_ch.protein_fasta, agt_ch.cds_fasta, vcf_ch)
*/
}

/*
workflow workflow_amr {
    take:
    contigs_ch
    
    main:
    //AMR
    //AMR1-ABRIcate
    abricate_ch = POST_ANALYSIS_ABRICATE(contigs_ch, params.organism)

    //AMR2-RESFINDER
    resfinder_ch = POST_ANALYSIS_AMRFINDER(contigs_ch)

}
*/

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