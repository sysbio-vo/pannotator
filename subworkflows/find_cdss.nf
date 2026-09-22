include { FIND_CDS } from '../modules/find_cds.nf'

workflow FIND_CDSS {
    take:
    assembly_ch // val(sample_id), path(assembly), path(bakta_db)

    main:
    cds_results = FIND_CDS(assembly_ch)
    // cds_results.view() // DEBUG

    emit:
    cds_results
}
