include { DETECT_PSEUDOGENES as DETECT_PSEUDOGENES_MODULE } from '../modules/detect_pseudogenes.nf'

workflow DETECT_PSEUDOGENES {
    take:
    manifest_file_and_bakta_db // tuple path(manifest_file), path(bakta_db)

    main:
    cds_with_pseudogenes = DETECT_PSEUDOGENES_MODULE(manifest_file_and_bakta_db)

    emit:
    cds_with_pseudogenes
}
