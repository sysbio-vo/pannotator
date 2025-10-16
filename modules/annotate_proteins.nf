process ANNOTATE_PROTEINS {
    tag "${proteins_fa.getBaseName()}"
    label "protein_annotation"
    publishDir "${params.outdir}", mode: 'symlink'

    input:
    path(proteins_fa)

    output:
    path("annotated_proteins_bakta/unique_proteins_annotation.faa")

    script:
    """
    bakta_proteins --db ${params.bakta_db} \
    --output annotated_proteins_bakta \
    --prefix unique_proteins_annotation \
    --threads 1 \
    ${proteins_fa}
    """
}