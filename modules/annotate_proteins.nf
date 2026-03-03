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