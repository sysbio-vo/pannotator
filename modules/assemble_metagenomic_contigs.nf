process ASSEMBLE_CONTIGS {
    tag "metagenomic"
    label "metagenomic"
    cpus { params.metaspades_threads }
    memory { params.metaspades_mem_cap }
    publishDir params.outdir, enabled: params.save_intermediate, mode: 'copy'

    input:
    tuple val(sample_id), path(sample_fwd), path(sample_rev) // paired-end reads assumbed by default

    output:
    tuple val(sample_id), path("megahit_out_${sample_id}/final.contigs.fa")

    script:
    // """
    // metaspades.py -1 ${sample_fwd} -2 ${sample_rev} \
    // -o metaspades_output \
    // -t ${params.metaspades_threads} \
    // -m ${params.metaspades_mem_cap}
    // """
    """
    megahit -1 ${sample_fwd} -2 ${sample_rev} \
    -o megahit_out_${sample_id} \
    -t ${params.metaspades_threads} \
    -m 17179869184
    """
}


process ASSEMBLE_CONTIGS_METASPADES {
    tag "metagenomic"
    label "metagenomic"
    cpus { params.metaspades_threads }
    memory { params.metaspades_mem_cap }
    publishDir params.outdir, enabled: params.save_intermediate, mode: 'copy'

    input:
    tuple val(sample_id), path(sample_fwd), path(sample_rev) // paired-end reads assumbed by default

    output:
    tuple val(sample_id), path("metaspades_out_${sample_id}/final.contigs.fa")

    script:
    // """
    // metaspades.py -1 ${sample_fwd} -2 ${sample_rev} \
    // -o metaspades_output \
    // -t ${params.metaspades_threads} \
    // -m ${params.metaspades_mem_cap}
    // """
    """
    megahit -1 ${sample_fwd} -2 ${sample_rev} \
    -o metaspades_out_${sample_id} \
    -t ${params.metaspades_threads} \
    -m 17179869184
    """
}