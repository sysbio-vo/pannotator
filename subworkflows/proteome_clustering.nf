include { CLUSTER_SEQS } from '../modules/cluster_sequences.nf'

workflow CLUSTER_PROTEOME {
    take:
    seqs_channel // path(seqs_file)

    main:
    CLUSTER_SEQS(seqs_channel)

    emit:
    CLUSTER_SEQS.out // all_seqs.fasta, cluster.tsv, rep_seq.fasta
}