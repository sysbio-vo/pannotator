process FIND_CDS {
    tag "cds_search"
    label "cds_search"
    label 'bakta'
    publishDir params.outdir, enabled: params.save_intermediate, mode: 'copy'

    input:
    tuple val(sample_id), path(assembly), path(bakta_db)

    output:
    tuple val(sample_id), path("CDSS_bakta/${sample_id}.cds-only.faa"), path("CDSS_bakta/${sample_id}.cds-only.pkl")

    script:
    meta = params.meta ? "--meta" : "" // Bakta (Pyrodigal) metagenome mode
    """
    bakta --db ${bakta_db} --cds-only  \
    ${meta} \
    --output CDSS_bakta  \
    --prefix ${sample_id} \
    --threads ${task.cpus} \
    ${assembly}
    """
}
