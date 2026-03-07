process FIND_CDS {
    tag "cds_search"
    label "cds_search"
    label 'bakta'
    publishDir params.outdir, enabled: params.save_intermediate, mode: 'copy'

    input:
    tuple path(assembly), path(bakta_db)

    output:
    tuple path("CDSS_bakta/${output_prefix}.cds-only.faa"), path("CDSS_bakta/${output_prefix}.cds-only.${params.serializer_ext}")

    script:
    output_prefix = assembly.getBaseName()
    """
    bakta --db ${bakta_db} --cds-only  \
    --output CDSS_bakta  \
    --prefix ${output_prefix} \
    --threads ${task.cpus} \
    ${assembly}
    """
}
