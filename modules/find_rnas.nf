process FIND_RNAS {
    tag { meta.tag }
    label "rnas_search"
    label 'bakta'
    publishDir params.outdir, enabled: params.save_intermediate, mode: 'copy'

    input:
    tuple val(meta), path(assemblies), path(bakta_db)

    output:
    tuple val(meta), path("${meta.tag}.rna-only.pkl")

    script:
    def individual_pickles = meta.asm_ids.collect { asm_id -> "RNAS_bakta/${asm_id}.rna-only.pkl" }.join(' ')
    """
    export TMPDIR=\$(mktemp -d)
    echo \$TMPDIR

    cleanup() {
        echo "Cleaning up \$TMPDIR..."
        rm -rf "\$TMPDIR"
    }
    trap cleanup EXIT

    # Loop through assemblies in batch running bakta on each
    for asm in ${assemblies}; do
        id=\$(basename "\$asm" | sed -E 's/(\\.[^.]+)+\$//')
        bakta --db ${bakta_db} --rna-only  \
          --output RNAS_bakta  \
          --prefix "\$id" \
          --threads ${task.cpus} \
        "\$asm"
    done

    # Concatenate to batch-level files
    merge_pickles.py \\
        --assembly_ids ${meta.asm_ids.join(',')} \\
        --out ${meta.tag}.rna-only.pkl \\
        ${individual_pickles}
    """
}
