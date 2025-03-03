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
include { ASSEMBLE                                            }     from '../bin/assemble/main'
include { PROKKA                                              }     from '../bin/anotations/prokka/main'
include { QUAST                                               }     from '../bin/qc/quast/main'
include { BUSCO                                               }     from '../bin/qc/busco/main'
include { MULTIQC_2 as POST_MULTIQC                           }     from '../bin/qc/multiqc/main_2' 
include { AMR as POST_ANALYSIS_ABRICATE                       }     from '../bin/AMR/abricate/main'
include { AMR_2 as POST_ANALYSIS_AMRFINDER                    }     from '../bin/AMR/AMRFinder/main'
include { ARIBA                                               }     from '../bin/mlst/main'



workflow assemble {
    preprocess_output = workflow_pre_process()
    amrprocess_output = workflow_amr( preprocess_output.contigs_ch, preprocess_output.fq_gz_reads_ch )
    postprocess_output = workflow_post_process( preprocess_output.busco_ch, preprocess_output.quast_all_ch )

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
   
    //Final Quality control after trimming
    fastq_ch_after = FASTQC_QUALITY_FINAL(trimmed_read_ch.trimmed_reads.map{it -> it[1]})

    //de novo assemble
    assemble_denovo_ch = ASSEMBLE(trimmed_read_ch.trimmed_reads)
    contigs_ch = assemble_denovo_ch.contigs
    scaffolds_ch = assemble_denovo_ch.scaffolds
    
    assemble_files_ch = contigs_ch
                .join(scaffolds_ch)
                
    quast_input_ch = assemble_files_ch.join(trimmed_read_ch.trimmed_reads)
    
    //PROKKA
    prokka_ch = PROKKA(contigs_ch)

    //BUSCO
    busco_ch = BUSCO(contigs_ch)

    //QUAST
    quast_ch = QUAST(quast_input_ch)
    quast_all_ch = quast_ch.direct_quast

    //MULTIQC
    multiqc_ch = MULTIQC(fastqc_ch_original.qc_zip.collect(), fastq_ch_after.qc_zip.collect())
    
    emit:
    contigs_ch
    fq_gz_reads_ch
    busco_ch
    quast_all_ch
}

workflow workflow_amr {
    take:
    contigs_ch
    fq_gz_reads_ch
    
    main:
    //AMR
    //AMR1-ABRIcate
    abricate_ch = POST_ANALYSIS_ABRICATE(contigs_ch)

    //AMR2-RESFINDER
    resfinder_ch = POST_ANALYSIS_AMRFINDER(contigs_ch)

    //MLST

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
    quast_all_ch

    main:
    multiqc_2_ch = POST_MULTIQC(quast_all_ch, busco_ch)

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