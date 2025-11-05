process MERGE_ANNOTATIONS {
    tag "merge_annotations"
    label "merge_annotations"

    publishDir params.outdir, mode: 'copy'

    input:
    path(cds_index)
    path(bulk_annotations)
    path(gff3_files)

    output:
    path('annotated_gff3')

    script:
    """
    map_annotations_to_samples.py \
        --gff3_dir . \
        --cds_index cds_index.json \
        --annotations bulk_protein_annotations.json \
        --out annotated_gff3
    """
}