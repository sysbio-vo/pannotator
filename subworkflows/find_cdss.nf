include { FIND_CDS } from '../modules/find_cds.nf'

workflow FIND_CDSS {
    take:
    indir // path(assembly)

    main:
    FIND_CDS(indir)
    
    emit:
    // outdir
    FIND_CDS.out
}
