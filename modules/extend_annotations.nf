process EXTEND_ANNOTATIONS {
    tag "extend_annotations"
    label "extend_annotations"

    publishDir params.outdir, mode: 'copy'

    input:
    path(clusters_pairs_tsv) // TODO: I would use tuple here
    path(all_seqs_fasta)
    path(bulk_annotations)

    output:
    path("bulk_protein_annotations_extended.json"), emit: bulk_annotations_extended

    script:
    """
    propagate_mmseqs_annotations.py \
    --tsv ${clusters_pairs_tsv} \
    --fasta ${all_seqs_fasta} \
    --json-in ${bulk_annotations} \
    --json-out bulk_protein_annotations_extended.json \
    --add-trace
    """
}
