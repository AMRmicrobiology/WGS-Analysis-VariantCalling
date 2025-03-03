#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Definir los parámetros con valores predeterminados
params.mode = params.mode ?: ''
params.input = params.input ?: ''
params.reference = params.reference ?: ''
params.outdir = params.outdir ?: 'results'

// Lista de modos válidos
def valid_modes = ['novo', 'reference', 'assemble']

// Convertir `--mode` en una lista (por si el usuario pasa varios workflows separados por comas)
def selected_modes = params.mode.split(',').collect { it.trim() }

// Verificar si todos los valores pasados en `--mode` son válidos
if( !selected_modes.every { it in valid_modes } ) {
    error "Invalid mode(s): '${params.mode}'. Please specify one or more of: 'novo', 'reference', 'assemble'."
}

// Incluir los sub-workflows desde la carpeta `subworkflow/`
include { novo } from './subworkflow/novo'
include { reference } from './subworkflow/reference'
include { assemble } from './subworkflow/assemble'

workflow {
    log.info """
    ==============================================
            WGS - N F   P I P E L I N E 
    ==============================================
    Running mode(s): ${selected_modes.join(', ')}
    Configuration environemnt:
    Out directory:             $params.outdir
    Fastq directory:           $params.input
    Reference directory:       $params.reference
    """

    selected_modes.each { mode ->
        switch (mode) {
            case 'novo':      novo(); break
            case 'reference': reference(); break
            case 'assemble':  assemble(); break
        }
    }
}