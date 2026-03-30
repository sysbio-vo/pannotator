process EXTEND_OR_GENERATE_AUXILIARY_DB {
    tag "auxiliary_db"
    label "auxiliary_db"

    publishDir file(params.auxiliary_db).parent, mode: 'copy'

    input:
    path(annotation_pkls)
    path(auxiliary_db) // auxiliary_db path, annotation pickles as collected list of paths

    output:
    path(auxiliary_db_name)

    script:
    auxiliary_db_name = auxiliary_db.name
    """
    generate_auxiliary_db.py \
        --annotation_pickles ${annotation_pkls.join(' ')} \
        --auxiliary_db ${auxiliary_db} \
        --updated_db_out ${auxiliary_db_name}
    """
}
