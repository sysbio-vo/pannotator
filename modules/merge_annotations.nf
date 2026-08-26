process MERGE_ANNOTATIONS {
    tag { meta.tag }
    label "merge_annotations"

    publishDir "${params.outdir}/predicted_cds_with_protein_function", enabled: params.save_intermediate, mode: 'copy'

    input:
    tuple val(meta), path(batch_cds_pkl)
    path(bulk_annotations)

    output:
    tuple val(meta), path("${out_pickle_dir}/*.cds-annotated.pkl"), emit: batch_annotated_pickles

    script:
    out_pickle_dir = "annotated_pkl"
    """
    merge_annotations_into_pkl.py \
        --pickle_in ${batch_cds_pkl} \
        --annotations ${bulk_annotations} \
        --pickle_out ${out_pickle_dir}
    """
}
