
process CLUSTER_SEQS {
    tag "clustering"
    label "clustering"

    input:
    tuple path(seqs_file) 

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