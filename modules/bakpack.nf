process BAKPACK_REDUCE {
    tag { meta.tag }
    label "bakpack_reduce"
    container "quay.io/sangerpathogens/bakpack:v0.5.0"

    publishDir "${params.outdir}/bakpack_reduced_annotations", mode: 'copy'

    input:
    tuple val(meta), path(json_tar), path(genomes_tar)

    output:
    path(out_fmt),                   emit: reduced_json_bakpack_bundle
    path("${meta.tag}.genomes.txt"), emit: assembly_list

    script:
    // no `def`: out_fmt must be visible to the output block
    if ( params.bakpack_reduce_outfmt == "zipped" ) {
        out_fmt = "${meta.tag}.reduced.json.tar.gz"
    } else if ( params.bakpack_reduce_outfmt == "json" ) {
        out_fmt = "${meta.tag}.reduced.json.dir"
    }
    """
    # one original assembly path per line
    printf '%s\n' ${meta.asm_paths.collect { "'${it}'" }.join(' ')} > ${meta.tag}.genomes.txt

    bakpack reduce --annotations ${json_tar} \
        --genomes ${genomes_tar} \
        -o ${out_fmt}
        
    rm ${genomes_tar} ${json_tar}
    """
}