process BUILD_INDEX {
    tag "Index-ReferenceGenome ${reference_id.name}"

    publishDir "${params.reference}/personal", mode: 'copy'

    input:
    path reference_id

    output:
    path "index_*_genome*", emit: index_data

    script:
    def index_prefix = "index_personal_genome"

    // Bowtie2    
    """
    bowtie2-build ${reference_id} ${index_prefix}
    """
}