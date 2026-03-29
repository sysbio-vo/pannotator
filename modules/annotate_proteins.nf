process ANNOTATE_PROTEINS {
    memory '10GB' // TODO: make versatile
    tag { proteins_fa.getBaseName() }
    label "protein_annotation"
    label 'bakta'
    publishDir params.outdir, enabled: params.save_intermediate, mode: 'copy'

    input:
    tuple path(proteins_fa), path(bakta_db)

    output:
    path("annotated_proteins_bakta/unique_proteins_annotation.json")

    script:
    """
    bakta_proteins --db ${bakta_db} \
    --output annotated_proteins_bakta \
    --prefix unique_proteins_annotation \
    --threads ${task.cpus} \
    ${proteins_fa}
    """
}

process ANNOTATE_WITH_AUX_DB {
    tag { proteins_fa.getBaseName() }
    label 'auxDB'
    publishDir params.outdir, enabled: params.save_intermediate, mode: 'copy'

    input:
    tuple path(proteins_fa), path(auxiliary_db)

    output:
    tuple path("auxiliary_annotations.json"), path("proteins_to_annotate.fa")

    script:
    """
    # perform annotation lookup in the auxiliary DB

    lookup_auxiliary_annotations.py \
    --auxiliary_db ${auxiliary_db} \
    --proteins_fa ${proteins_fa} \
    --out auxiliary_annotations.json
    --remaining_proteins_filename proteins_to_annotate.fa
    """
}

process MERGE_BULK_ANNOTATIONS {
    tag "auxDB"
    publishDir params.outdir, enabled: params.save_intermediate, mode: 'copy'

    input:
    path(annotation_jsons)

    output:
    path('merged_bulk_protein_annotation.json')

    script:
    """
    merge_bulk_json_annotations.py \
    -i ${annotation_jsons.join(' ')}
    --out merged_bulk_protein_annotation.json
    """
}

process PARSE_BAKTA_JSON_ANNOTATIONS {
    tag { bakta_proteins_json_annot.getBaseName() }
    publishDir params.outdir, enabled: params.save_intermediate, mode: 'copy'

    input:
    path(bakta_proteins_json_annot)

    output:
    path("bulk_protein_annotations.json")

    script:
    """
    parse_annotation_json.py -a ${bakta_proteins_json_annot} \
    -o bulk_protein_annotations.json
    """
}