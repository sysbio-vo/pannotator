include { BUILD_COORDS_INDEX } from '../modules/build_coords_index.nf'

workflow BUILD_COORDS_INDEX_WF {
    take:
    cds_dir_path
    outdir

    main:
    BUILD_COORDS_INDEX(cds_dir_path, outdir)

    emit:
    BUILD_COORDS_INDEX.out
}
