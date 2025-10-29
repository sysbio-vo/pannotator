process BUILD_COORDS_INDEX {
    tag "build_cds_index"
    label "build_cds_index"

    publishDir params.outdir, mode: 'copy'

    input:
    path(cds_dir)

    output:
    path('cds_index.json') 

    script:
    """
    build_coords_index.py ${cds_dir} cds_index.json
    """
}
