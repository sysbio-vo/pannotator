include { FIND_CDS } from '../modules/find_cds.nf'

workflow FIND_CDSS {
    take:
    assemblies_channel // path(assembly)

    main:
    FIND_CDS(assemblies_channel)
    FIND_CDS.out
        .set { found_cdss }

    emit:
    found_cdss
}