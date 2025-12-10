process BAKTA {
    tag "BAKTA annotation for ${sample_id}"
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        "docker://${params.bakta.docker}" :
        params.bakta.docker }"

    publishDir "${params.outdir}/2-Assembly/Annotations", mode: 'copy', saveAs: { filename ->
        if (filename.endsWith(".gff3")) {
            return "bakta/${sample_id}/${sample_id}.gff3"
        } else if (filename.endsWith(".faa")) {
            return "bakta/${sample_id}/${sample_id}.faa"
        } else if (filename.endsWith(".fna")) {
            return "bakta/${sample_id}/${sample_id}.fna"
        } else if (filename.endsWith(".gbff")) {
            return "bakta/${sample_id}/${sample_id}.gbff"
        } else if (filename.endsWith(".txt")) {
            return "bakta/${sample_id}/${sample_id}.txt"
        } else if (filename.endsWith(".json")) {
            return "bakta/${sample_id}/${sample_id}.json"
        } else {
            return null
        }
    }

    input:
    tuple val(sample_id), path(assembly_file)

    output:
    path "annotations_${sample_id}/${sample_id}.gff3", emit: bakta_gff3
    path "annotations_${sample_id}/${sample_id}.faa", emit: bakta_faa
    path "annotations_${sample_id}/${sample_id}.ffn", emit: bakta_ffn
    path "annotations_${sample_id}/${sample_id}.gbff", emit: bakta_gbff
    path "annotations_${sample_id}/${sample_id}.txt", emit: bakta_txt
    path "annotations_${sample_id}/${sample_id}.json", emit: bakta_json
    tuple val(sample_id), path("annotations_${sample_id}/${sample_id}.gff3"), path("annotations_${sample_id}/${sample_id}.fna"), emit: conv_gff

    
    script:

    """
    amrfinder_update --force_update --database /data/db-light/amrfinderplus-db

    bakta --db /data/db-light --threads ${task.cpus} --keep-contig-headers --output annotations_${sample_id} ${assembly_file}
    """
}