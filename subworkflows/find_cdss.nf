include { FIND_CDS } from '../modules/find_cds.nf'

workflow FIND_CDSS {
    take:
    indir // path(assembly)
    outdir 

    main:
    FIND_CDS(indir, outdir)
    
    emit:
    // outdir
    FIND_CDS.out
}
