nextflow.enable.dsl = 2

checkInputParams()

reference = file("${params.reference}")

log.info """\

WGS - N F   P I P E L I N E 
            FOR 
    C L I N I C A L   A N A L Y S I S
==============================================
Configuration environment:
    Out directory:             $params.outdir
    Fastq directory:           $params.input
    Reference directory:       $params.reference
"""
    .stripIndent()

// Incluir los sub-workflows según el modo seleccionado
if (params.mode == 'novo') {
    include { novo } from './subworkflow/novo' 
} else if (params.mode == 'reference') {
    include { reference } from './subworkflow/reference'
} else {
    error "Invalid mode: ${params.mode}. Please specify 'novo' or 'reference'."
}

// Definir el workflow principal
workflow {
    if (params.mode == 'novo') {
        novo()  // Llamar al workflow 'novo' que ha sido incluido
    } else if (params.mode == 'reference') {
        reference()  // Llamar al workflow 'reference' que ha sido incluido
    }
}


////////////////////////////////////////////////////////////////////////////////
// FUNCTIONS                                                                  //
////////////////////////////////////////////////////////////////////////////////

def checkInputParams() {
    // Check required parameters and display error messages
    boolean fatal_error = false

    if (!params.input) {
        log.warn("You need to provide a fastqDir (--fastqDir) or a bamDir (--bamDir)")
        fatal_error = true
    }

    // Check required parameters based on the mode selected
    if (params.mode == 'novo') {
        // In clinical mode, no user-supplied reference is required
        log.info("Using wildtype reference for clinical workflow")
    } else if (params.mode == 'reference') {
        // In reference mode, a user-supplied reference is required
        if (!params.reference) {
            log.warn("You need to provide a genome reference (--reference) for the reference workflow")
            fatal_error = true
        }
    } else {
        log.warn("Invalid mode: ${params.mode}. Please specify 'clinical' or 'reference'.")
        fatal_error = true
    }
    
    if (fatal_error) {
        error("Required parameters are missing.")
    }
}
