process MERGE_FASTA {
    tag "aggregate"
    label "aggregate"
    publishDir params.outdir, mode: 'symlink'

    input:
    path(aa_seqs)

    output:
    path('all_proteins.faa')

    script:
    """
    cat ./*faa > all_proteins.faa
    """
}


process CLUSTER_SEQS {
    tag "clustering"
    label "clustering"
    publishDir params.outdir, mode: 'symlink'

    input:
    path(seqs_file)

    output:
    tuple path("${params.mmseqs_clusterPrefix}_all_seqs.fasta"), path("${params.mmseqs_clusterPrefix}_cluster.tsv"), path("${params.mmseqs_clusterPrefix}_rep_seq.fasta")

    script:
    """
    mmseqs ${params.mmseqs_command} \
    ${seqs_file} \
    ${params.mmseqs_clusterPrefix} \
    ${params.mmseqs_tmpDir} \
    ${params.mmseqs_args} 
    """
}