include { ANNOTATE_PROTEINS as ANNOTATE_PROTEINS_MODULE; 
          ANNOTATE_PROTEINS_WITH_AUXILIARY_DB; 
          PARSE_BAKTA_JSON_ANNOTATIONS } from '../modules/annotate_proteins.nf'

workflow ANNOTATE_PROTEINS {
    take:
    unique_proteins_channel // path(proteins_fa), path(bakta_db), path(auxiliary_db)

    main:
    if ( params.auxiliary_db ) {
        ANNOTATE_PROTEINS_WITH_AUXILIARY_DB(unique_proteins_channel)
        ANNOTATE_PROTEINS_WITH_AUXILIARY_DB.out
            .set { protein_annotations }
    } else {
        ANNOTATE_PROTEINS_MODULE(unique_proteins_channel)
        ANNOTATE_PROTEINS_MODULE.out
            .set { protein_annotations }
    }

    PARSE_BAKTA_JSON_ANNOTATIONS(protein_annotations)
    PARSE_BAKTA_JSON_ANNOTATIONS.out
        .set { bulk_annotations }

    emit:
    protein_annotations
    bulk_annotations
}