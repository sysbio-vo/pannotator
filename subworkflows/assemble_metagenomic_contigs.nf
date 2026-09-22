include { ASSEMBLE_CONTIGS } from '../modules/assemble_metagenomic_contigs.nf'

workflow ASSEMBLE_MAGS {
    take:
    indir // tuple val(sample_id), path(assembly)

    main:
    reads_str = "${params.indir}/*{${params.sample_1_suffix},${params.sample_2_suffix}}*${params.infile_extension}"
    println "${reads_str} DEBUG reads_str" 

    read_pairs_ch = Channel.fromFilePairs(reads_str, checkIfExists: true)
    read_pairs_ch.map { meta, reads -> [meta, reads[0], reads[1]] }
        .set {read_pairs_ch}
    // read_pairs_ch.view() // DEBUG

    contigs = ASSEMBLE_CONTIGS(read_pairs_ch)

    emit:
    contigs
}