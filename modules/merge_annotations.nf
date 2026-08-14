process MERGE_ANNOTATIONS {
    tag "merge_annotations"
    label "merge_annotations"

    publishDir "${params.outdir}/protein_annotations", enabled: params.save_intermediate, mode: 'copy'

    input:
    path(batch_cds_pkl)
    path(bulk_annotations)

    output:
    path("annotated_pkl/*.cds-annotated.pkl"), emit: batch_annotated_pickles

    script:
    """
    merge_annotations_into_pkl.py \
        --pickle_in . \
        --annotations ${bulk_annotations} \
        --pickle_out annotated_pkl
    """
}
