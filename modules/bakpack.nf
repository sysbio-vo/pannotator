process BAKPACK_REDUCE {
    tag { meta.tag }
    label "bakpack_reduce"
    container "quay.io/sangerpathogens/bakpack:v0.5.0"

    publishDir (
        "${params.outdir}/bakpack_reduced_annotations",
        pattern: "*.reduced.json.tar.gz",
        mode: 'copy'
    )
    publishDir (
        "${params.outdir}/bakpack_reduced_annotations",
        pattern: "*.genomes.tar.gz",
        mode: 'copy'
    )

    input:
    tuple val(meta), path(json_tar), path(genomes_tar)

    output:
    path("${meta.tag}.reduced.json.tar.gz"), emit: reduced_json_bundle
    path("${meta.tag}.genomes.tar.gz"),      emit: genome_bundle

    script:
    """
    bakpack reduce --annotations ${json_tar} \
        --genomes ${genomes_tar} \
        -o ${meta.tag}.reduced.json.tar.gz
    
    tar -czf ${meta.tag}.genomes.tar.gz ${genomes_tar} 
    """
}