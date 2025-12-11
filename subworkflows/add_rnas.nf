include { FIND_RNAS } from '../modules/find_rnas.nf'

workflow ADD_RNAS {
    take:
    indir 

    main:
    annotated_rna = FIND_RNAS(indir)  

    emit:
    annotated_rna
}
