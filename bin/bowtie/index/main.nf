process BUILD_INDEX {
    tag "Index-ReferenceGenome ${reference_id.name}"

    publishDir "${params.reference}/personal/index", mode: 'copy'

    input:
    tuple val(sample_id), path(reference_id)

    output:
    path("index_personal_genome.*.bt2"), emit: index_files


    script:
    def index_prefix = "index_personal_genome"

    // Bowtie2    
    """
    bowtie2-build ${reference_id} ${index_prefix}
    """
}