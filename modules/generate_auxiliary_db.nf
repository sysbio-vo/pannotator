process EXTEND_OR_GENERATE_AUXILIARY_DB {
    tag "auxiliary_db"
    label "auxiliary_db"
    label "bakta" // signifies utilisation of Bakta functions in the executable

    publishDir params.auxiliary_db ? file(params.auxiliary_db).parent : "", mode: 'copy'

    input:
    path(annotation_pkls)
    path(cds_annotation_pkls)
    tuple path(auxiliary_db), path(bulk_annotation_before_filtering) // auxiliary_db path, annotation pickles as collected list of paths

    output:
    path(auxiliary_db_name)

    script:
    auxiliary_db_name = auxiliary_db.name
    """
    generate_auxiliary_db.py \
        --cds_annotation_pickles ${cds_annotation_pkls.join(' ')} \
        --sorf_annotation_pickles ${annotation_pkls.join(' ')} \
        --auxiliary_db ${auxiliary_db} \
        --bulk_annotation_before_filtering ${bulk_annotation_before_filtering} \
        --updated_db_out ${auxiliary_db_name}
    """
}
