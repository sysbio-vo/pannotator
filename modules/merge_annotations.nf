process MERGE_ANNOTATIONS {
    tag "merge_annotations"
    label "merge_annotations"

    publishDir params.outdir, enabled: params.save_intermediate, mode: 'copy'

    input:
    path(cds_pkl_files)
    path(bulk_annotations)

    output:
    path("annotated_pkl/*.pkl"), emit: annotated_pickles
    path("annotated_pkl"), emit: annotated_dir

    script:
    """
    merge_annotations_into_pkl.py \
        --pickle_folder . \
        --annotations ${bulk_annotations} \
        --out annotated_pkl
    """
}
