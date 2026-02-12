include { CLUSTER_SEQS } from '../modules/cluster_sequences.nf'

workflow CLUSTER_PROTEOME {
    take:
    seqs_channel // path(aa_seqs)

    main:

    seqs_channel
        .filter( ~/.*\.faa/ )
        .collectFile( name: 'all_proteins.faa', newLine: true)
        .set { all_proteins_ch }
            
    CLUSTER_SEQS(all_proteins_ch)

    emit:
    CLUSTER_SEQS.out // all_seqs.fasta, cluster.tsv, rep_seq.fasta
}