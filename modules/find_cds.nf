process FIND_CDS {
    tag { meta.tag }
    label "cds_search"
    label 'bakta'
    publishDir params.outdir, enabled: params.save_intermediate, mode: 'copy'

    input:
    tuple val(meta), path(assemblies), path(bakta_db)

    output:
    tuple val (meta), path("${meta.tag}.cds-only.faa"), path("${meta.tag}.cds-only.pkl")

    script:
    """
    # Loop through assemblies in batch running bakta on each
    for asm in ${assemblies}; do
        id=\$(basename "\$asm" | sed -E 's/(\\.[^.]+)+\$//')
        bakta --db ${bakta_db} --cds-only  \
          --output CDSS_bakta  \
          --prefix "\$id" \
          --threads ${task.cpus} \
        "\$asm"
    done
    
    # Concatenate fastas
    cat CDSS_bakta/*.faa > ${meta.batch_id}.cds-only.faa

    # TODO: Need to do this for pickles too, but in a script maybe?
    # cat CDSS_bakta/*.pkl > ${meta.batch_id}.cds-only.pkl
    """
}
