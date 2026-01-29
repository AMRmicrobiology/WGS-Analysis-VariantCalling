process PROKKA {
    tag "PROKKA ANNOTATION"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "docker://${params.prokka.docker}" :
        params.prokka.docker }"

    publishDir "${params.outdir}/2-Assembly/3-Annotations", mode: 'copy', saveAs: { filename ->
        if (filename.endsWith(".gff")) {
            return "prokka/${sample_id}/${sample_id}.gff"
        } else if (filename.endsWith(".faa")) {
            return "prokka/${sample_id}/${sample_id}.faa"
        } else if (filename.endsWith(".fna")) {
            return "prokka/${sample_id}/${sample_id}.fna"
        } else {
            return null
        }
    }

    input:
    tuple val(sample_id), path(assembly_file)

    output:
    tuple val(sample_id), path("annotations_${sample_id}/${sample_id}_wildtype.gff"), emit: prokka_tuple_gff
    path "annotations_${sample_id}/${sample_id}_wildtype.gff", emit: prokka_gff
    path "annotations_${sample_id}/${sample_id}_wildtype.faa", emit: prokka_faa
    path "annotations_${sample_id}/${sample_id}_wildtype.fna", emit: prokka_fna

    script:
    """
    prokka --outdir annotations_${sample_id} --prefix ${sample_id}_wildtype --kingdom Bacteria --centre AMRmicrobiology --compliant ${assembly_file}
    """
}