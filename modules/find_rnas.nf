process FIND_RNAS {
    tag "rnas_search"
    label "rnas_search"
    label 'bakta'
    publishDir params.outdir, enabled: params.save_intermediate, mode: 'copy'

    input:
    tuple path(assembly), path(bakta_db)

    output:
    path("RNAS_bakta/${output_prefix}.rna-only.pkl")

    script:
    output_prefix = assembly.getBaseName()
    """
    bakta --db ${bakta_db} --rna-only  \
    --output RNAS_bakta  \
    --prefix ${output_prefix} \
    --threads ${task.cpus} \
    ${assembly}
    """
}
