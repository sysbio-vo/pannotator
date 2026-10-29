process BAKPACK_BUILD {
    tag { meta.tag }
    container "quay.io/sangerpathogens/bakpack:v0.4.0"

    publishDir (
        "${params.outdir}/bakpack_annotations",
        pattern: "*.bakpack",
        mode: 'copy'
    )
    publishDir (
        "${params.outdir}/bakpack_annotations",
        pattern: "*.tar.xz",
        mode: 'copy'
    )

    input:
    tuple val(meta), path(tar_json_annotations), path(genomes)

    output:
    tuple val(meta), path("${meta.tag}_annotations.bakpack"), emit: bakpack
    tuple val(meta), path("${meta.tag}_annotations.tar.xz}"), path("${meta.tag}_genomes.tar.xz"), emit: tar_bundles

    script:

    """
    bakpack build \
        --annotations ${tar_json_annotations} \
        --genomes ${genomes} \
        --xz-threads ${task.cpus} \
        --output ${meta.tag}_annotations.bakpack
    """
}