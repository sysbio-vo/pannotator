include { FIND_CDS } from '../modules/find_cds.nf'
include { COLLECT_CDS_FILES } from '../modules/collect_cds.nf'

workflow FIND_CDSS {
    take:
    indir // path(assembly)

    main:
    cds_results = FIND_CDS(indir)
    collected_files = cds_results
        .flatten()
        .collect()
    
    COLLECT_CDS_FILES(collected_files)

    emit:
    COLLECT_CDS_FILES.out
}
