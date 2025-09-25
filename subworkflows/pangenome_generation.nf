include { BUILD_PANGENOME } from '../modules/build_pangenome.nf' 
include { CLUSTER_PROTEOME } from 'proteome_clustering.nf'

workflow GENERATE_PANGENOME {
    take:
    assemblies_channel // path(assembly)

    main:
    FIND_CDSS(assemblies_channel) | CLUSTER_PROTEOME | BUILD_PANGENOME
        .set { pangenome_index }

    emit:
    pangenome_index
}