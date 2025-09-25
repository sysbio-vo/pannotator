workflow FIND_CDS {
    tag "cds_search"
    label "cds_search"

    input:
    path(assembly)

    output:
    path(found_cdss)

    script:
    """
    python find_cds_bakta.py \
    -i ${assembly} \
    -o ${found_cdss}
    """
}