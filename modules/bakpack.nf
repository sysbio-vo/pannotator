process BAKPACK_REDUCE {
    tag { meta.tag }
    label "bakpack_reduce"
    container "quay.io/sangerpathogens/bakpack:v0.4.0"

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
    # untar the archives
    tar -xf ${json_tar}
    tar -xf ${genomes_tar}

    # run bakpack reduce on all of the jsons
    for json in *.json; 
        do bakpack reduce \${json} \${json/.json/.fa} -o \${json/.json/_reduced.json};
    done

    # tar up the reduced files
    tar -czf ${meta.tag}.reduced.json.tar.gz *_reduced.json && \
    tar -czf ${meta.tag}.genomes.tar.gz *.fa && \

    rm *.json *.fa
    """
}