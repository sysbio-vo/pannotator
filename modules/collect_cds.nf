process COLLECT_CDS_FILES {
    tag "collect_cds_files"

    input: 
    path(cds_files)

    output:
    path("collected_cds_files")

    script:
    """
    mkdir -p collected_cds_files
    cp ${cds_files} collected_cds_files/
    """
}