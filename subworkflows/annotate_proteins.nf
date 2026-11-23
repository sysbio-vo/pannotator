include { ANNOTATE_PROTEINS as ANNOTATE_PROTEINS_MODULE;
          ANNOTATE_WITH_AUX_DB as ANNOTATE_WITH_AUX_DB_MODULE;
          ANNOTATE_WITH_AUX_DB_URL;
          MERGE_BULK_ANNOTATIONS;
          PARSE_BAKTA_JSON_ANNOTATIONS } from '../modules/annotate_proteins.nf'

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


workflow ANNOTATE_WITH_AUX_DB {
    take:
    unique_proteins_channel // path(proteins_fa), path(bakta_db)
    hmm_ch // path(user_hmms)
    prots_ch // path(user_proteins)

    main:
    def auxIsUrl = params.auxiliary_db.toString().startsWith('http://') || params.auxiliary_db.toString().startsWith('https://')

    unique_proteins_channel
        .multiMap { proteins_fa, bakta_db ->
            proteins_fa: proteins_fa
            bakta_db: bakta_db
        }
        .set { unique_proteins }

    if (auxIsUrl) {
        ANNOTATE_WITH_AUX_DB_URL(
            unique_proteins.proteins_fa.map { proteins_fa -> tuple(proteins_fa, params.auxiliary_db.toString()) }
        )
        aux_lookup = ANNOTATE_WITH_AUX_DB_URL.out
    } else {
        ANNOTATE_WITH_AUX_DB_MODULE(
            unique_proteins.proteins_fa.combine(Channel.of(file(params.auxiliary_db)))
        )
        aux_lookup = ANNOTATE_WITH_AUX_DB_MODULE.out
    }

    aux_lookup
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

    ANNOTATE_PROTEINS(remaining_proteins_to_annotate, hmm_ch, prots_ch)

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