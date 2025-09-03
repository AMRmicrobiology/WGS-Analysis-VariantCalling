nextflow.enable.dsl=2

checkInputParams()

reference         = file("${params.reference}")


//Call all the sub-work
include { FASTQC_QUALITY as FASTQC_QUALITY_ORIGINAL           }     from '../bin/qc/fastqc/main'
include { TRIMMING                                            }     from '../bin/trimming/main'
include { FASTQC_QUALITY as FASTQC_QUALITY_FINAL              }     from '../bin/qc/fastqc/main'
include { PREPARE_KRAKEN_DB                                   }     from '../bin/kraken/prepare_db'
include { KRAKEN;SEQTK_PRUNE                                  }     from '../bin/kraken/main'
include { MULTIQC                                             }     from '../bin/qc/multiqc/main'
include { BAKTA                                               }     from '../bin/anotations/bakta/main'
include { EXTRACT_CDS_FROM_BAKTA			                  }     from '../bin/anotations/bakta/main_2'
include { BUILD_INDEX_1                                       }     from '../bin/bowtie/index/main_bwa'
include { BUILD_INDEX as PERSONAL_GENOME_INDEX                }     from '../bin/bowtie/index/main'
include { PERSONAL_GENOME_MAPPING                             }     from '../bin/bowtie/mapping/main'
include { MARKDUPLICATE                                       }     from '../bin/gatk/picard/markduplicate/main'
include { ADDORREPLACE                                        }     from '../bin/gatk/picard/addorreplace/main'
include { HAPLOTYPECALLER                                     }     from '../bin/gatk/haplotype/main_1'
include { GENOTYPE as GENOTYPE_ANALYSIS                       }     from '../bin/gatk/genotype/main'
include { ALIGN as NORMALISE_DATA                             }     from '../bin/gatk/Filter/align'
include { FILTER_VARIANTS as FILTER_VARIANTS_PARAM            }     from '../bin/gatk/Filter/main'
include { DECOMPRESS_VCF                                      }     from '../bin/snpeff/main_2'
include { SNPEFF			                                  }     from '../bin/snpeff/main_3'

workflow reference {
    krakenprocess_output = workflow_kraken_process()
    preprocess_output = workflow_pre_process(krakenprocess_output.DB_CH)
    
    postprocess_output = workflow_post_process(preprocess_output.fq_gz_reads_ch,
        preprocess_output.prune_reads_ch)
    
}

workflow workflow_kraken_process {
    db_ready_ch = PREPARE_KRAKEN_DB()
    DB_CH= db_ready_ch.db_ready

    emit:
    DB_CH
}

workflow workflow_pre_process {
    take:
    DB_CH

    main:
    // Quality control and Index build
    read_ch = Channel.fromFilePairs(params.input, size: 2)
    
    fastqc_ch_original= FASTQC_QUALITY_ORIGINAL(read_ch.map{it -> it[1]})
    
    // Trimming process
    trimmed_read_ch = TRIMMING(read_ch)
    fq_gz_reads_ch = trimmed_read_ch.trimmed_reads

    //KRAKEN
    READS_DB_CH = fq_gz_reads_ch.combine(DB_CH)
                .map { sample_id, reads_pair, db_dir ->
                def (r1, r2) = reads_pair
                tuple (sample_id, [r1, r2], db_dir)
    }

    kraken_ch = KRAKEN (READS_DB_CH)

    //Final Quality control after trimming
    fastq_ch_after = FASTQC_QUALITY_FINAL(trimmed_read_ch.trimmed_reads.map{it -> it[1]})
    
    //PRUNNING
    fastq_prunning_ch = fq_gz_reads_ch.join(kraken_ch.keep_ids).map {
        sample_id,reads_pair, keep_ids ->
        def (r1, r2) = reads_pair
        tuple (sample_id, [r1, r2], keep_ids)
    }
    
    prune_ch = SEQTK_PRUNE(fastq_prunning_ch)
    prune_reads_ch = prune_ch.pruned_reads

    //MULTIQC
    multiqc_ch = MULTIQC(fastqc_ch_original.qc_zip.collect(), fastq_ch_after.qc_zip.collect())



    emit:
    fq_gz_reads_ch
    prune_reads_ch

}

workflow workflow_post_process {

    take:
    prune_reads_ch
    fq_gz_reads_ch

    main:

    //Reference Genome INDEX
    personal_ref_ch = Channel.fromPath(params.personal_ref)
    reference_ch = personal_ref_ch.map {
        ref -> 
        def sample_id = file(ref).baseName
        def ref_file = file(ref)
        return tuple (sample_id, ref_file)
    }

    personal_index_ch = PERSONAL_GENOME_INDEX(reference_ch)
    
    //mapping process- Mapping used Specie ref. genome, include samtools sorted
    mapping_input_ch = prune_reads_ch.combine(personal_index_ch)

    specie_mapping_ch = PERSONAL_GENOME_MAPPING(mapping_input_ch)

    //Add groups and add or replace group
    bam_ch = specie_mapping_ch.map {
        tupla -> 
        def sample_id = tupla [0]
        def bam_path = tupla [2]
        return tuple (sample_id, bam_path)
    }

    gatk_mark_ch = ADDORREPLACE(bam_ch)
    
    //Marckduplicate
    gatk_add_ch = MARKDUPLICATE(gatk_mark_ch)
    
    //HAPLOTYPECALLER realignment consistently

    haplotype_ch = gatk_add_ch.dedup_bam.combine(reference_ch.map {id_reference, reference -> tuple(id_reference, reference) })
    .set { all_samples_ch }

    gatk_haplotype_ch= HAPLOTYPECALLER(all_samples_ch) 
    
    //GenotypeCaller
    gatk_genotype_ch = GENOTYPE_ANALYSIS (gatk_haplotype_ch.out_files.combine(reference_ch.map {id_reference, reference -> tuple(id_reference, reference) }))
    
    //Align
    //This tool takes a VCF file, left-aligns the indels and trims common bases from indels, leaving them with a minimum representation.
    //The same indel can often be placed at multiple positions and still represent the same haplotype.
    //We are going to take the optionally splits multiallelic sites into biallelics and left-aligns individual alleles.
    aligns_and_normalized_ch = NORMALISE_DATA (gatk_genotype_ch.combine(reference_ch.map {id_reference, reference -> tuple(id_reference, reference) }))
    
    //VatiantFilter
    //Filter the VCF using the parametres to get a hight quality and cover in SNPs and INDELS "QUAL || MQ || DP ".
    //all the parametres could be changen it, depends of the data.
    variant_filter_ch = FILTER_VARIANTS_PARAM (aligns_and_normalized_ch.combine(reference_ch.map {id_reference, reference -> tuple(id_reference, reference) }))
   
    // Decompress VCF
    vcf_ch = DECOMPRESS_VCF(variant_filter_ch.compl_vcf)

    // BAKTA PROCESS BUILD A GFF AND CDS OF REFERENCE
    // PARAMS GFF PROVIDE OR NOT FORM THE CUSTOMER
    // Select GFF source (BAKTA or custom)
    if (params.custom_gff3 && params.custom_gff3 != 'null') {
        log.info "Using custom GFF3 file provided by the user: ${params.custom_gff3}"
        gff3_ch = Channel.value(file(params.custom_gff3))
    } else {
        log.info "No custom GFF3 file provided — running BAKTA to generate it from the reference"
        gff_first_ch = BAKTA(reference_ch)
        gff3_ch = gff_first_ch.bakta_gff3
    }
    
    // Combine channels for SNPEFF
    vcf_gff_combined_ch = vcf_ch.combine(gff3_ch)
    vcf_gff_ref_combined_ch = vcf_gff_combined_ch.combine(reference_ch)

    snpeff_input_ch = vcf_gff_ref_combined_ch.map { entry ->
        def (sample_id, vcf_path, gff3_path, ref_id, ref_fasta) = entry
        return tuple(
            gff3_path,
            ref_id,
            ref_fasta,
            params.genome_name_db,
            sample_id,
            vcf_path
        )
    }

    snpeff_ch = SNPEFF(snpeff_input_ch)
      
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