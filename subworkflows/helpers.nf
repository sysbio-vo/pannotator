

workflow BATCHER {
    take:
    infiles // path(assembly_file)
    bakta_db // path(bakta_db)
    buffer_size // val(int)

    main:
    infiles
    .map { asm -> tuple([id: sampleIdFromName(asm.name)], asm) }
    .buffer (size: buffer_size, remainder: true)
    .toList()
    .flatMap { batches -> 
        batches.withIndex().collect {batch, idx ->
            def meta = [
                batch_id: idx,
                asm_ids: batch.collect { asm_tuple -> asm_tuple[0].id },
                tag: batch.size() == 1 ? batch [0][0].id : "batch_${idx}" // adapts to per asm / per batch of assemblies
            ] // TODO: check if remainder is 1, if that batch is tagged differently to the rest
            tuple(meta, batch.collect { asm_tuple -> asm_tuple[1] })
        }
    }
    .combine(bakta_db)
    .set { batches_and_bakta_db }

    emit:
    batches_and_bakta_db
}