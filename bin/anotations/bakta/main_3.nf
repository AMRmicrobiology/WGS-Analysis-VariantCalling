process BAKTA {
    tag "BAKTA annotation for ${sample_code}" 

    container {
        workflow.containerEngine == 'singularity' ?
            "docker://${params.bakta_db.docker}" :
            params.bakta_db.docker
    }
    publishDir "${params.outdir}/2-Assembly/3-Annotations", mode: 'copy', saveAs: { filename ->
        // Tu lógica de publicación original
        if (filename.endsWith(".gff3")) return "bakta/${sample_code}/${sample_code}.gff3"
        else if (filename.endsWith(".faa")) return "bakta/${sample_code}/${sample_code}.faa"
        else if (filename.endsWith(".fna")) return "bakta/${sample_code}/${sample_code}.fna"
        else if (filename.endsWith(".gbff")) return "bakta/${sample_code}/${sample_code}.gbff"
        else if (filename.endsWith(".txt")) return "bakta/${sample_code}/${sample_code}.txt"
        else if (filename.endsWith(".json")) return "bakta/${sample_code}/${sample_code}.json"
        else if (filename.endsWith(".ffn")) return "bakta/${sample_code}/${sample_code}.ffn" // Aseguramos el FFN
        else return null
    }

    input:
    tuple val(sample_code), path(assembly_file)
    path db_directory

    output:
    tuple val(sample_code), path("annotations_${sample_code}/${sample_code}_consensus_wrapped.gff3"), emit: bakta_gff3
    path "annotations_${sample_code}/${sample_code}_consensus_wrapped.faa", emit: bakta_faa
    path "annotations_${sample_code}/${sample_code}_consensus_wrapped.ffn", emit: bakta_ffn
    path "annotations_${sample_code}/${sample_code}_consensus_wrapped.gbff", emit: bakta_gbff
    path "annotations_${sample_code}/${sample_code}_consensus_wrapped.txt", emit: bakta_txt
    path "annotations_${sample_code}/${sample_code}_consensus_wrapped.json", emit: bakta_json

    script:
    """
    OUTPUT_DIR="annotations_${sample_code}"
    
    bakta --db ${db_directory}/db-light \
          --threads ${task.cpus} \
          --keep-contig-headers \
          --skip-plot \
          --prefix ${sample_code}_consensus_wrapped \
          --output \${OUTPUT_DIR} \
          ${assembly_file}

    """
}
