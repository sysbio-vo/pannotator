process DETECT_PSEUDOGENES {
    tag "cds_pseudogenes"
    label "cds_psedogenes"
    label 'bakta'
    publishDir params.outdir, enabled: params.save_intermediate, mode: 'copy'

    input:
    tuple path(manifest_file), path(bakta_db)

    output:
    path("cds_with_pseudogenes/*with_pseudogenes.pkl")

    script:
    // output_prefix = assembly.getBaseName()
    // bakta_db_arg = params.bakta_db ? "--db ${params.bakta_db}" : ""
    """
    bakta_pseudo_bulk \
    --db ${bakta_db} \
    --output cds_with_pseudogenes \
    --prefix pseudogenes \
    ${manifest_file}
    """
}