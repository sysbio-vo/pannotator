include { ANNOTATE_PROTEINS as ANNOTATE_PROTEINS_MODULE } from '../modules/annotate_proteins.nf'

workflow ANNOTATE_PROTEINS {
    take:
    unique_proteins_channel // path(proteins_fa)

    main:
    ANNOTATE_PROTEINS_MODULE(unique_proteins_channel)
    ANNOTATE_PROTEINS_MODULE.out
        .set { protein_annotations }

    emit:
    protein_annotations
}