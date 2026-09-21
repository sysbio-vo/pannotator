include { FIND_CDS } from '../modules/find_cds.nf'

workflow FIND_CDSS {
    take:
    batches_and_bakta_db // tuple(meta, list_of_assembly_paths, bakta_db_path)

    main:
    batched_cds_results = FIND_CDS(batches_and_bakta_db)

    emit:
    batched_cds_results
}
