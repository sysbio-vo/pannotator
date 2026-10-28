process BAKPACK_BUILD {


    publishDir (
        "${params.outdir}/bakpack_annotations",
        pattern: "${meta.tag}_annotations.bakpack",
        mode: 'copy'
    )
    
    container "quay.io/sangerpathogens/bakpack:v0.4.0"

    input:
    tuple val(meta), path(tar_json_annotations), path(genomes)

    output:
    tuple val(meta)

    script:

    """
    bakpack build \
        --annotations ${tar_json_annotations} \
        --genomes ${genomes} \
        --xz-threads ${task.cpus} \
        --output ${meta.tag}_annotations.bakpack
    """
}