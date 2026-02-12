process SORF_EXTRA {
    tag "${sample_id} sorf_extra_search"
    label "sorf_extra_search"
    publishDir params.outdir, enabled: params.save_intermediate, mode: 'copy'

    input:
    tuple val(sample_id), path(assembly), path(cds_pkl), path(rna_pkl)

    output:
    tuple val(sample_id), path("${sample_id}")

    script:
    """
    bakta --db ${params.bakta_db} --sorf-extra  \
        --cds-pickle ${cds_pkl} \
        --rna-pickle ${rna_pkl} \
        --output ${sample_id}  \
        --prefix ${sample_id} \
        --threads ${task.cpus} \
        ${assembly}
    """
}