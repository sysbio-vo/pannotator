process BAKPACK_REDUCE {
    tag { meta.tag }
    label "bakpack_reduce"
    container "quay.io/sangerpathogens/bakpack:v0.5.0"

    publishDir (
        "${params.outdir}/bakpack_reduced_annotations",
        pattern: "*.reduced.json.bpk",
        mode: 'copy'
    )
    // publishDir (
    //     "${params.outdir}/bakpack_reduced_annotations",
    //     pattern: "*.genomes.tar.gz",
    //     mode: 'copy'
    // )
    publishDir (
        "${params.outdir}/bakpack_reduced_annotations",
        pattern: "*.genomes.txt",
        mode: 'copy'
    )

    input:
    tuple val(meta), path(json_tar), path(genomes_tar)

    output:
    path("${meta.tag}.reduced.json.bpk"),    emit: reduced_json_bakpack_bundle
    path("${meta.tag}.genomes.txt"),      emit: assembly_list

    script:
    """
    # one original assembly path per line
    printf '%s\n' ${meta.asm_paths.collect { "'${it}'" }.join(' ')} > ${meta.tag}.genomes.txt


    bakpack reduce --annotations ${json_tar} \
        --genomes ${genomes_tar} \
        -o ${meta.tag}.reduced.json.bpk
        
    rm ${genomes_tar} ${json_tar}
    """
}