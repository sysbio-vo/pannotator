process BAKPACK_BUILD {
    container "quay.io/sangerpathogens/bakpack:v0.4.0"

    publishDir (
        "${params.outdir}/bakpack_annotations",
        pattern: "${meta.tag}_annotations.bakpack",
        mode: 'copy'
    )
    publishDir (
        "${params.outdir}/bakpack_annotations",
        pattern: "${meta.tag}.*.tar.xz",
        mode: 'copy'
    )

    input:
    tuple val(meta), path(tar_json_annotations), path(genomes)

    output:
    tuple val(meta), "${meta.tag}_annotations.bakpack"

    script:

    """
    bakpack build \
        --annotations ${tar_json_annotations} \
        --genomes ${genomes} \
        --xz-threads ${task.cpus} \
        --output ${meta.tag}_annotations.bakpack
    """
}