process FIND_RNAS {
    tag "rnas_search"
    label "rnas_search"
    label 'bakta'
    publishDir params.outdir, enabled: params.save_intermediate, mode: 'copy'

    input:
    tuple val(sample_id), path(assembly), path(bakta_db)

    output:
    tuple val(sample_id), path("RNAS_bakta/${sample_id}.rna-only.pkl")

    script:
    """
    export TMPDIR=\$(mktemp -d)
    echo \$TMPDIR

    cleanup() {
        echo "Cleaning up \$TMPDIR..."
        rm -rf "\$TMPDIR"
    }
    trap cleanup EXIT

    bakta --db ${bakta_db} --rna-only  \
    --verbose \
    --output RNAS_bakta  \
    --prefix ${sample_id} \
    --threads ${task.cpus} \
    ${assembly}
    """
}
