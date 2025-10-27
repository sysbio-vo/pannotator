process BUILD_COORDS_INDEX {
    tag "build_cds_index"
    label "build_cds_index"

    publishDir { outdir }, mode: 'copy'

    input:
    val(cds_dir)
    val(outdir) 

    output:
    path('cds_index.json') 

    script:
    """
    build_coords_index.py ${cds_dir} cds_index.json
    """
}
