process DETECT_PSEUDOGENES {
    tag { meta.tag }
    tag "cds_pseudogenes"
    label "cds_pseudogenes"
    label 'bakta'
    publishDir "${params.outdir}/pseudogenes", enabled: params.save_intermediate, mode: 'copy'

    input:
    tuple path(batch_pickle), path(bakta_db)

    output:
    path("cds_with_pseudogenes/*with_pseudogenes.pkl")

    script:
    // output_prefix = assembly.getBaseName()
    // bakta_db_arg = params.bakta_db ? "--db ${params.bakta_db}" : ""
    """
    # Loop through assemblies in batch running bakta on each
    for asm in ${assemblies}; do
        id=\$(basename "\$asm" | sed -E 's/(\\.[^.]+)+\$//')
        bakta_pseudo_bulk \
        --db ${bakta_db} \
        --output cds_with_pseudogenes \
        --prefix pseudogenes \
        --batch_pickle ${batch_pickle}
    done
    """
}