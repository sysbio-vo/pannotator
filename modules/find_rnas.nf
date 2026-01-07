process FIND_RNAS {
    tag "rnas_search"
    label "rnas_search"
    publishDir params.outdir, enabled: params.save_intermediate, mode: 'copy'

    input:
    path(assembly)

    output:
    path("RNAS_bakta/${output_prefix}.rna-only.pkl")

    script:
    output_prefix = assembly.getBaseName()
    """
    bakta --db ${params.bakta_db} --rna-only  \
    --output RNAS_bakta  \
    --prefix ${output_prefix} \
    --threads ${task.cpus} \
    ${assembly}
    """
}