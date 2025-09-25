process BUILD_PANGENOME {
    tag "pangenome_construction"
    label "pangenome_construction"

    input:
    path(clustered_seqs)

    output:
    path(pangenome_index)

    script:
    """
    python build_pangenome.py \
    -i ${clustered_seqs} \
    -o ${pangenome_index}
    """
}