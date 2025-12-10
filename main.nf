/*
  ============================================================
 __        ______ ____        _                _           _        
 \ \      / / ___/ ___|      / \   _ __   __ _| |_   _ ___(_)___    
  \ \ /\ / / |  _\___ \     / _ \ | '_ \ / _` | | | | / __| / __|   
   \ V  V /| |_| |___) |   / ___ \| | | | (_| | | |_| \__ \ \__ \   
    \_/\_/  \____|____/   /_/   \_\_| |_|\__,_|_|\__, |___/_|___/   
 __     __         _             _      ____     |___/_             
 \ \   / /_ _ _ __(_) __ _ _ __ | |_   / ___|__ _| | (_)_ __   __ _ 
  \ \ / / _` | '__| |/ _` | '_ \| __| | |   / _` | | | | '_ \ / _` |
   \ V / (_| | |  | | (_| | | | | |_  | |__| (_| | | | | | | | (_| |
    \_/ \__,_|_|  |_|\__,_|_| |_|\__|  \____\__,_|_|_|_|_| |_|\__, |
                                                              |___/
      N F   P I P E L I N E - WGS_ANALYSIS_VARIANT_CALLING

  Illumina sequencing data WGS/Variant Calling Pipeline - Nextflow
  ============================================================

  Author:        Jimmy Lucas and Roger de Pedro Jové
  Description:   Nextflow pipeline for whole-genome sequencing (WGS)
                 analysis and variant calling in bacterial genomes 
                 using Illumina data, supporting de novo assembly and 
                 reference-based analysis.
  Version:       1.0.0

  ============================================================
*/

nextflow.enable.dsl = 2

if (params.help) {
    printHelp()
    exit 0
}

checkInputParams()

reference = file("${params.reference}")

log.info """
 __        ______ ____        _                _           _        
 \\ \\      / / ___/ ___|      / \\   _ __   __ _| |_   _ ___(_)___    
  \\ \\ /\\ / / |  _\\___ \\     / _ \\ | '_ \\ / _` | | | | / __| / __|   
   \\ V  V /| |_| |___) |   / ___ \\| | | | (_| | | |_| \\__ \\ \\__ \\   
    \\_/\\_/  \\____|____/   /_/   \\_\\_| |_|\\__,_|_|\\__, |___/_|___/   
 __     __         _             _      ____     |___/_             
 \\ \\   / /_ _ _ __(_) __ _ _ __ | |_   / ___|__ _| | (_)_ __   __ _ 
  \\ \\ / / _` | '__| |/ _` | '_ \\| __| | |   / _` | | | | '_ \\ / _` |
   \\ V / (_| | |  | | (_| | | | | |_  | |__| (_| | | | | | | | (_| |
    \\_/ \\__,_|_|  |_|\\__,_|_| |_|\\__|  \\____\\__,_|_|_|_|_| |_|\\__, |
                                                              |___/ 

==============================================
N F   P I P E L I N E - WGS_ANALYSIS_VARIANT  
==============================================
Configuration environment:
    Pipeline mode:             $params.mode
    Fastq directory:           $params.input
    Profile:                   $workflow.profile

"""
    .stripIndent()

// Subworkflows 

if (params.mode == 'assemble') {
    include { assemble } from "$projectDir/subworkflow/assemble" 
} else if (params.mode == 'reference') {
    include { reference } from "$projectDir/subworkflow/reference"
} else if (params.mode == 'novo') {
    include { novo} from "$projectDir/subworkflow/novo"
} else {
    error "Invalid mode: ${params.mode}. Please specify 'assemble' ,'reference' or 'novo'"
}

// Definir el workflow principal
workflow {
    if (params.mode == 'assemble') {
        assemble()  
    } else if (params.mode == 'reference') {
        reference()
    } else if (params.mode == 'novo') {
        novo()
    }
}

////////////////////////////////////////////////////////////////////////////////
// FUNCTIONS                                                                  //
////////////////////////////////////////////////////////////////////////////////

def printHelp() {
    def readmeFile = file("${projectDir}/README.md")
    def printSection = false

    if (readmeFile.exists()) {
        log.info "\n"
        readmeFile.eachLine { line ->
            // Start printing when we hit the Usage header
            if (line.contains("Usage: nextflow run main.nf [--help] [--mode VAR] [--input VAR] [--genome_name_db VAR] [--wildtype_code VAR] [--outdir VAR] [--personal_ref VAR] [--custom_gff3 VAR] [--organism VAR] [--cut_front VAR] [--cut_tail VAR] [--cut_mean_quality VAR] [--length_required VAR] [--mrsa] -[-qual_snp VAR] [--qual_indel VAR] [--bakta_db_define VAR] [--db_select VAR] [--abricate_db VAR] [-w VAR] [-profile VAR]")) {
                printSection = true
            }
            // Stop printing when we hit the next major header (Output)
            if (line.contains("## Output")) {
                printSection = false
            }
            
            // Print the line if we are inside the section
            if (printSection) {
                log.info line
            }
        }
        log.info "\n"
    } else {
        log.warn "README.md not found in ${projectDir}"
    }
}

def checkInputParams() {
    // Check required parameters and display error messages
    boolean fatal_error = false

    if (!params.input) {
        log.warn("You need to provide a valid input directory with --input")
        fatal_error = true
    }
    if (!params.mode) {
        log.warn("You need to provide a valid mode with --mode (assemble, novo, reference)")
        fatal_error = true
    }
    if( params.mode == 'reference' && !params.personal_ref ) {
        log.warn "You need to provide a valid personal reference with --personal_ref when using reference mode"
        fatal_error = true
    }
    if( !['docker','singularity','conda'].contains( workflow.profile ) ) {
        log.warn "You need to provide a valid profile with -profile (docker, singularity, conda)"
        fatal_error = true
    }
    if( params.mode == 'assemble' && !params.organism ) {
        log.warn "You need to provide a valid organism with --organism when using assemble mode"
        fatal_error = true
    }
    if (fatal_error) {
        error "Missing one or more required parameters"
    }
    
}