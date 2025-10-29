include { MERGE_FASTA; CLUSTER_SEQS } from '../modules/cluster_sequences.nf'

workflow CLUSTER_PROTEOME {
    take:
    seqs_dir_channel // path(aa_seqs_dir)

    main:
    MERGE_FASTA(seqs_dir_channel)
    MERGE_FASTA.out
        .set { all_proteins_ch }
    CLUSTER_SEQS(all_proteins_ch)

    emit:
    CLUSTER_SEQS.out // all_seqs.fasta, cluster.tsv, rep_seq.fasta
}