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

process ANNOTATE_PROTEINS_WITH_AUXILIARY_DB {
    tag { proteins_fa.getBaseName() }
    label "protein_annotation"
    label 'bakta'
    publishDir params.outdir, enabled: params.save_intermediate, mode: 'copy'

    input:
    tuple path(proteins_fa), path(bakta_db), path(auxiliary_db)

    output:
    path("annotated_proteins_bakta/unique_proteins_annotation.json")

    script:
    """
    # perform annotation lookup in the auxiliary DB first (create a new one if doesn't exist)
    lookup_auxiliary_annotations.py \
    --auxiliary_db ${auxiliary_db} \
    --proteins_fa ${proteins_fa} \
    --out auxiliary_annotations.json
    --remaining_proteins_filename proteins_to_annotate.fa


    # annotate proteins that were not found in the auxiliary DB
    bakta_proteins --db ${bakta_db} \
    --output annotated_proteins_bakta \
    --prefix bakta_proteins_annotation \
    --threads ${task.cpus} \
    proteins_to_annotate.fa


    # merge JSON annotations
    merge_anno_jsons.py \
    --anno_jsons annotated_proteins_bakta/bakta_proteins_annotation.json auxiliary_annotations.json \
    --out annotated_proteins_bakta/unique_proteins_annotation.json
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