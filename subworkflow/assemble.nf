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
include { KRAKEN;SEQTK_PRUNE                                  }     from '../bin/kraken/main'
include { ASSEMBLE                                            }     from '../bin/assemble/main'
include { FILTER_CONTIGS                                      }     from '../bin/qc/polish/filter'
include { ALIGMENT_PILON;PILON_POLISH                         }     from '../bin/qc/polish/main'
include { PROKKA                                              }     from '../bin/anotations/prokka/main'
include { QUAST                                               }     from '../bin/qc/quast/main'
include { BUSCO                                               }     from '../bin/qc/busco/main'
include { MULTIQC_2 as POST_MULTIQC                           }     from '../bin/qc/multiqc/main_2' 
include { MRSA                                                }     from '../bin/mrsa/main'
include { SCCMEC                                              }     from '../bin/mrsa/main'
include { AMR as POST_ANALYSIS_ABRICATE                       }     from '../bin/AMR/abricate/main'
include { AMR_2 as POST_ANALYSIS_AMRFINDER                    }     from '../bin/AMR/AMRFinder/main'
include { ARIBA                                               }     from '../bin/mlst/main'
include { MLST                                                }     from '../bin/mlst/main_2'


workflow assemble {
    preprocess_output = workflow_pre_process()
    amrprocess_output = workflow_amr( preprocess_output.accurance_fasta_ch, preprocess_output.fq_gz_reads_ch )
    postprocess_output = workflow_post_process( preprocess_output.busco_ch, preprocess_output.quast_ch )
    if (params.mrsa) {
        mrsaprocess_output = workflow_mrsa(preprocess_output.accurance_fasta_ch)
    }
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
    prokka_ch = PROKKA(accurance_fasta_ch)

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
    
    emit:
    accurance_fasta_ch
    fq_gz_reads_ch
    busco_ch
    quast_ch
}

workflow workflow_amr {
    take:
    accurance_fasta_ch
    fq_gz_reads_ch
    
    main:
    //AMR
    //AMR1-ABRIcate
    abricate_ch = POST_ANALYSIS_ABRICATE(accurance_fasta_ch)

    //AMR2-RESFINDER
    resfinder_ch = POST_ANALYSIS_AMRFINDER(accurance_fasta_ch)

    //MLST FAST RAW DATA- ARIBA

    def organism_schemes_ch = Channel.fromPath('organisms_list.txt')
        .splitText()
        .map { line -> line.trim() }
        .filter { it.startsWith(params.organism) }
        .map { scheme -> tuple(params.organism, scheme) }
        .unique()

    def combined_ch = fq_gz_reads_ch.combine(organism_schemes_ch)

    ariba_ch = ARIBA(combined_ch)

}

workflow workflow_post_process {

    take:
    busco_ch
    quast_ch

    main:
    multiqc_2_ch = POST_MULTIQC(quast_ch, busco_ch)

}

workflow workflow_mrsa {
    take:
    accurance_fasta_ch

    main:
    
    //MRSA

    mrsa_ch = MRSA (accurance_fasta_ch)
    sccmec_ch = SCCMEC(accurance_fasta_ch)

    //MLST
    MLST(accurance_fasta_ch)
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