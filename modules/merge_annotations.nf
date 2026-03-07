process MERGE_ANNOTATIONS {
    tag "merge_annotations"
    label "merge_annotations"

    publishDir params.outdir, mode: 'copy'

    input:
    path(cds_pkl_files)
    path(bulk_annotations)

    output:
    path("annotated_pkl/*.${params.serializer_ext}"), emit: annotated_pickles
    path("annotated_pkl"), emit: annotated_dir

    script:
    """
    merge_annotations_into_pkl.py \
        --infile_ext ${params.serializer_ext} \
        --pickle_folder . \
        --annotations bulk_protein_annotations.json \
        --out annotated_pkl
    """
}
