include { ANNOTATE_PROTEINS as ANNOTATE_PROTEINS_MODULE;
          ANNOTATE_WITH_AUX_DB as ANNOTATE_WITH_AUX_DB_MODULE;
          MERGE_BULK_ANNOTATIONS;
          PARSE_BAKTA_JSON_ANNOTATIONS } from '../modules/annotate_proteins.nf'

workflow ANNOTATE_PROTEINS {
    take:
    unique_proteins_channel // path(proteins_fa), path(bakta_db)

    main:
    ANNOTATE_PROTEINS_MODULE(unique_proteins_channel)
    ANNOTATE_PROTEINS_MODULE.out
        .set { protein_annotations }

    PARSE_BAKTA_JSON_ANNOTATIONS(protein_annotations)
    PARSE_BAKTA_JSON_ANNOTATIONS.out
        .set { bulk_annotations }

    emit:
    protein_annotations
    bulk_annotations
}


workflow ANNOTATE_WITH_AUX_DB {
    take:
    unique_proteins_channel // path(proteins_fa), path(bakta_db), path(auxiliary_db)

    main:
    unique_proteins_channel
        .multiMap { proteins_fa, bakta_db, auxiliary_db -> 
            bulk_proteins_and_aux_db: tuple(proteins_fa, auxiliary_db)
            bakta_db: bakta_db
        }
        .set { unique_proteins }

    ANNOTATE_WITH_AUX_DB_MODULE(unique_proteins.bulk_proteins_and_aux_db)

    ANNOTATE_WITH_AUX_DB_MODULE.out
        .multiMap { aux_anno, remaining_proteins ->
            bulk_aux_anno: aux_anno
            proteins_to_annotate: remaining_proteins
        }
        .set { aux_anno_and_remaining_proteins }
    
    aux_anno_and_remaining_proteins.proteins_to_annotate
        .map { fasta -> tuple(fasta, fasta.countFasta()) }
        .set { proteins_to_annotate_with_count }


    proteins_to_annotate_with_count.filter { file, protein_count -> protein_count > 0 }
        .map { file, protein_count -> file }
        .set { filtered_proteins_to_annotate }

    filtered_proteins_to_annotate
        .combine(unique_proteins.bakta_db)
        .set { remaining_proteins_to_annotate }

    ANNOTATE_PROTEINS(remaining_proteins_to_annotate)

    ANNOTATE_PROTEINS
        .out
        .bulk_annotations
        .set { novel_annotated_proteins_json }

    novel_annotated_proteins_json
        .concat(aux_anno_and_remaining_proteins.bulk_aux_anno)
        .collect()
        .set { all_bulk_annotation_jsons }

    MERGE_BULK_ANNOTATIONS(all_bulk_annotation_jsons)

    MERGE_BULK_ANNOTATIONS
        .out
        .set { bulk_annotations }

    emit:
    bulk_annotations
}