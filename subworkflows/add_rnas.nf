include { FIND_RNAS } from '../modules/find_rnas.nf'

workflow ADD_RNAS {
    take:
    indir // tuple(meta, list_of_assembly_paths, bakta_db_path)

    main:
    annotated_rna = FIND_RNAS(batches_and_bakta_db)  

    emit:
    annotated_rna
}
