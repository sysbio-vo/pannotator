include { FIND_CDS } from '../modules/find_cds.nf'

workflow FIND_CDSS {
    take:
    indir // path(assembly)

    main:
    cds_results = FIND_CDS(indir)
    collected_files = cds_results
        .flatten()
        .collect()    

    emit:
    collected_files
}
