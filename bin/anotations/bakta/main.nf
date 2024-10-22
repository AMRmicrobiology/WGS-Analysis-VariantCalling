process BAKTA {
    tag "BAKTA ANNOTATIONS"
    container params.bakta.docker

    publishDir "${params.outdir}", mode: 'copy', saveAs: { filename ->
        if (filename.endsWith(".gff3")) {
            return "7-bakta/${sample_id}/${sample_id}.gff3"
        } else if (filename.endsWith(".faa")) {
            return "7-bakta/${sample_id}/${sample_id}.faa"
        } else if (filename.endsWith(".fna")) {
            return "7-bakta/${sample_id}/${sample_id}.fna"
        } else {
            return null
        }
    }

    input:
    tuple val(sample_id), path(assembly_file)

    output:
    path "annotations_${sample_id}/${sample_id}.gff3", emit: bakta_gff3
    path "annotations_${sample_id}/${sample_id}.faa", emit: bakta_faa
    path "annotations_${sample_id}/${sample_id}.fna", emit: bakta_fna

    script:

    """
    amrfinder_update --force_update --database /data/db-light/amrfinderplus-db

    bakta --db /data/db-light --threads ${task.cpus} --keep-contig-headers --output annotations_${sample_id} ${assembly_file}

    """
}