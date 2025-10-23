process ANNOTATE_PROTEINS {
    memory '10GB'
    tag "${proteins_fa.getBaseName()}"
    label "protein_annotation"
    publishDir "${params.outdir}", mode: 'symlink'

    input:
    path(proteins_fa)

    output:
    path("annotated_proteins_bakta/unique_proteins_annotation.json")

    script:
    """
    bakta_proteins --db ${params.bakta_db} \
    --output annotated_proteins_bakta \
    --prefix unique_proteins_annotation \
    --threads 1 \
    ${proteins_fa}
    """
}

process PARSE_BAKTA_JSON_ANNOTATIONS {
    publishDir "${params.outdir}", mode: 'symlink'

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