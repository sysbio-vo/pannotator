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

    main:

    // predict pseudogenes using annotated pickle objects
    
    annotated_samples
        .multiMap { meta, pkl ->
            meta_ch: meta
            pkl_ch: pkl
        }
        .set { split_samples }

    // TODO: Nextflow caching doesn't work well with this approach
    // if a single new sample is added, this whole subworkflow reruns
    split_samples.pkl_ch
        .flatten()
        .buffer (size: params.batch_of_batches_size, remainder: true)
        .map { it.join('\n') } // TODO: is this crutch REALLY neccessary to collect paths to files in a txt files instead of their contents?
        .collectFile( name: 'annotated_cds_manifest.txt', newLine: true, sort: true )
        .set { manifest_file }

    // bakta_pseudo_bulk processes every original batch's pickles in one bulk
    // call, so there's no single per-batch meta left to carry through; meta
    // is only used for the process `tag {}` label, so any one of them does
    split_samples.meta_ch
        .first()
        .set { bulk_meta }

    bulk_meta
        .combine(manifest_file)
        .combine(bakta_db)
        .set { manifest_file_and_bakta_db }

    DETECT_PSEUDOGENES(manifest_file_and_bakta_db)

    annotated_samples_updated = DETECT_PSEUDOGENES.out
        .flatten()


    emit:
    annotated_samples_updated
}