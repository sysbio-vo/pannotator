process FIND_CDS {
    tag { meta.tag }
    label "cds_search"
    label 'bakta'
    publishDir "${params.outdir}/predicted_cds", enabled: params.save_intermediate, mode: 'copy'

    input:
    tuple val(meta), path(assemblies), path(bakta_db)

    output:
    tuple val (meta), path("${meta.tag}.cds-only.faa"), path("${meta.tag}.cds-only.pkl")

    script:
    def individual_pickles = meta.asm_ids.collect { asm_id -> "CDSs_bakta/${asm_id}.cds-only.pkl" }.join(' ')
    """
    # Loop through assemblies in batch running bakta on each
    for asm in ${assemblies}; do
        id=\$(basename "\$asm" | sed -E 's/(\\.[^.]+)+\$//')
        bakta --db ${bakta_db} --cds-only  \
          --output CDSs_bakta  \
          --prefix "\${id}" \
          --threads ${task.cpus} \
          --force \
        "\$asm"
    done
    
    # Concatenate outputs to batch-level files
    cat CDSs_bakta/*.faa > ${meta.tag}.cds-only.faa
    manage_pkls.py batch \\
        --sample-ids ${meta.asm_ids.join(',')} \\
        --output ${meta.tag}.cds-only.pkl \\
        ${individual_pickles} \
    # remove the intermediate directory to save space
    rm -rf CDSs_bakta/
    """
}
