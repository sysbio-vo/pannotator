include { DETECT_PSEUDOGENES as DETECT_PSEUDOGENES_MODULE } from '../modules/detect_pseudogenes.nf'

workflow DETECT_PSEUDOGENES {
    take:
    manifest_file_and_bakta_db // tuple path(manifest_file), path(bakta_db)

    main:
    cds_with_pseudogenes = DETECT_PSEUDOGENES_MODULE(manifest_file_and_bakta_db)

    emit:
    cds_with_pseudogenes
}

workflow DETECT_PSEUDOGENES_OPTIONAL {
    take:
    annotated_samples
    bakta_db
    bakta_db_type

    main:

    if ( bakta_db_type == 'full' ) {
        // predict pseudogenes using annotated pickle objects

        // TODO: Nextflow caching doesn't work well with this approach
        // if a single new sample is added, this whole subworkflow reruns
        annotated_samples
            .flatten()
            .map { it -> "${it}" } // TODO: is this crutch REALLY neccessary to collect paths to files in a txt files instead of their contents? 
            .collectFile( name: 'annotated_cds_manifest.txt', newLine: true, sort: true )
            .set { manifest_file }
        
        manifest_file
            .combine(bakta_db)
            .set { manifest_file_and_bakta_db }

        DETECT_PSEUDOGENES(manifest_file_and_bakta_db)

        annotated_samples_updated = DETECT_PSEUDOGENES.out
            .flatten()
    } else {
        annotated_samples_updated = annotated_samples
            .flatten()
    }


    emit:
    annotated_samples_updated
}