process FIND_CDS {
    tag "cds_search"
    label "cds_search"
    publishDir params.outdir, mode: 'symlink'

    container 'quay.io/d_goryslavets/bakta_pannotator:1.11.4-pannotator.1.2'

    input:
    path(assembly)

    output:
    tuple path("CDSS_bakta/${output_prefix}.cds-only.faa"), path("CDSS_bakta/${output_prefix}.cds-only.gff3")

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
