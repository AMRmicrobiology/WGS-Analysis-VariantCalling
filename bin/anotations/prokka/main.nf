process PROKKA {
    tag "PROKKA ANNOTATION"

    publishDir "${params.outdir}", mode: 'copy', saveAs: { filename ->
        if (filename.endsWith(".gff")) {
            return "6-prokka/${sample_id}/${sample_id}.gff"
        } else if (filename.endsWith(".faa")) {
            return "6-prokka/${sample_id}/${sample_id}.faa"
        } else if (filename.endsWith(".fna")) {
            return "6-prokka/${sample_id}/${sample_id}.fna"
        } else {
            return null
        }
    }

    input:
    tuple val (sample_id), path(assembly_file)

    output:
    path "annotations_${sample_id}/${sample_id}_wildtype.gff", emit: prokka_gff
    path "annotations_${sample_id}/${sample_id}_wildtype.faa", emit: prokka_faa
    path "annotations_${sample_id}/${sample_id}_wildtype.fna", emit: prokka_fna

    script:
    """
    prokka --outdir annotations_${sample_id} --prefix ${sample_id}_wildtype --kingdom Bacteria ${assembly_file}
    """
}