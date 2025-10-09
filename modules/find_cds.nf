process FIND_CDS {
    tag "cds_search"
    label "cds_search"
    publishDir "${params.outdir}", mode: 'symlink'

    input:
    path(assembly)

    output:
    path("CDSS_bakta/${output_prefix}.cds-only.faa")

    script:
    output_prefix = assembly.getBaseName()
    // bakta_db_arg = params.bakta_db ? "--db ${params.bakta_db}" : ""
    """
    bakta --db ${params.bakta_db} --cds-only  \
    --output CDSS_bakta  \
    --prefix ${output_prefix} \
    --threads 1 \
    ${assembly}
    """
}