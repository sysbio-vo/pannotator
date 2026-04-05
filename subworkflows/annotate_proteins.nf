include { ANNOTATE_PROTEINS as ANNOTATE_PROTEINS_MODULE; PARSE_BAKTA_JSON_ANNOTATIONS } from '../modules/annotate_proteins.nf'

workflow ANNOTATE_PROTEINS {
    take:
    unique_proteins_channel // path(proteins_fa), path(bakta_db)
    hmm_ch //path(user_hmms)
    prots_ch // path(user_proteins)

    main:
    ANNOTATE_PROTEINS_MODULE(unique_proteins_channel, hmm_ch, prots_ch)
    ANNOTATE_PROTEINS_MODULE.out
        .set { protein_annotations }

    PARSE_BAKTA_JSON_ANNOTATIONS(protein_annotations)
    PARSE_BAKTA_JSON_ANNOTATIONS.out
        .set { bulk_annotations }

    emit:
    protein_annotations
    bulk_annotations
}