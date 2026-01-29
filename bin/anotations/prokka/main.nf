process PROKKA {
    tag "PROKKA annotation for ${sample_id}"

        container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "docker://${params.prokka.docker}" :
        params.prokka.docker }"

    publishDir "${params.outdir}/2-assemble/3-annotations", mode: 'copy', saveAs: { filename ->
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
    path "annotations_${sample_id}/${sample_id}.gff", emit: prokka_gff
    tuple val(sample_id), path "annotations_${sample_id}/${sample_id}.gff", emit: prokka_tuple_gff
    path "annotations_${sample_id}/${sample_id}.faa", emit: prokka_faa
    path "annotations_${sample_id}/${sample_id}.fna", emit: prokka_fna

    script:
    """
    prokka --outdir annotations_${sample_id} --prefix ${sample_id} --kingdom Bacteria --compliant ${assembly_file}
    
    """
}