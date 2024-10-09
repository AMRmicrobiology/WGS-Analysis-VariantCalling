process PROKKA {
    tag "PROKKA ANOTATION"

    input:
    tuple val(sample_id), path(assembly_file)

    output:
    path "annotations/${sample_id}_wildtype.gff", emit: prokka_file

    script:

    """
    mkdir -p annotations
    prokka --outdir annotations --prefix ${sample_id}_wildtype --kingdom Bacteria ${assembly_file}
    """
}