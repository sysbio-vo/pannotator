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

    // TODO: Nextflow caching doesn't work well with this approach
    // if a single new sample is added, this whole subworkflow reruns
    annotated_samples
        .buffer (size: params.batch_of_batches_size, remainder: true)
        .toList()
        .flatMap { chunks -> chunks.withIndex(1).collect { chunk, bob_num -> tuple("bob_${bob_num}", chunk) } }
        .set { bob_chunks } // (bob_tag, [[meta, pkl], ...])

    bob_chunks
        .map { bob_tag, chunk -> tuple(bob_tag, chunk.collect { meta, pkl -> pkl }.flatten().join('\n')) } // TODO: is this crutch REALLY neccessary to collect paths to files in a txt files instead of their contents?
        .collectFile( newLine: true ) { bob_tag, content -> [ "${bob_tag}_annotated_cds_manifest.txt", content ] }
        .map { manifest_file -> tuple(manifest_file.baseName - '_annotated_cds_manifest', manifest_file) }
        .set { bob_manifest_files }

    // keep the original per-batch metas so they survive the batch-of-batches regrouping
    bob_chunks
        .map { bob_tag, chunk -> tuple(bob_tag, chunk.collect { meta, pkl -> meta }) }
        .set { bob_metas } // (bob_tag, [original batch meta, ...])

    bob_manifest_files
        .join(bob_metas)
        .map { bob_tag, manifest_file, batch_metas -> tuple([tag: bob_tag, batch_metas: batch_metas], manifest_file) }
        .combine(bakta_db)
        .set { manifest_file_and_bakta_db }

    DETECT_PSEUDOGENES(manifest_file_and_bakta_db)

    // After batches of batches match original batches with original meta maps
    annotated_samples_updated = DETECT_PSEUDOGENES.out
        .flatMap { bob_meta, pkls ->
            def pkl_list = pkls instanceof List ? pkls : [pkls]
            pkl_list.collect { pkl ->
                def original_meta = bob_meta.batch_metas.find { bm -> pkl.name.startsWith("${bm.tag}.") }
                assert original_meta != null : "Could not match pseudogene output ${pkl.name} back to an original batch meta in ${bob_meta.tag}"
                tuple(original_meta, pkl)
            }
        }

    emit:
    annotated_samples_updated
}