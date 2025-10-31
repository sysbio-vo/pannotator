process MERGE_ANNOTATIONS {
    tag "merge_annotations"
    label "merge_annotations"

    publishDir params.outdir, mode: 'copy'

    input:
    path(cds_index)
    path(bulk_annotations)

    output:
    path('merged_annotations.json')

    script:
    """
    merge_annotations.py ${cds_index} ${bulk_annotations} -o merged_annotations.json
    """

}