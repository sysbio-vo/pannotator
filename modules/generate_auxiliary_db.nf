process EXTEND_OR_GENERATE_AUXILIARY_DB {
    tag "auxiliary_db"
    label "auxiliary_db"

    publishDir params.outdir, enabled: params.save_intermediate, mode: 'copy'

    input:
    tuple path(auxiliary_db), path(annotation_gff3s) // auxiliary_db.json, annotation_gff3s list of paths

    output:
    path("auxiliary_db.json"), emit: auxiliary_db

    script:
    """
    generate_alignment_db_from_result_gff3s.py \
        --auxiliary_db ${auxiliary_db}
        --annotations ${annotation_gff3s} \
        --out auxiliary_db.json
    """
}
