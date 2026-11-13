process SORF_EXTRA {
    tag { meta.tag }
    label "sorf_extra_search"
    label 'bakta'

    publishDir (
        "${params.outdir}/final_annotations",
        pattern: params.bundle_gff3 ? "*.gff3.tar.gz" : "*.gff3",
        mode: 'copy'
    )
    publishDir (
        "${params.outdir}/final_annotations",
        pattern: "SORFs_bakta/*.json",
        mode: 'copy',
        enabled: { !params.bundle_gff3 }
    )
    publishDir (
        "${params.outdir}/final_annotations",
        pattern: "*.sorf-extra.pkl",
        mode: 'copy',
        enabled: { params.save_intermediate }
    )

    input:
    tuple val(meta), path(assemblies), path(cds_pkl), path(rna_pkl), path(bakta_db)

    output:
    tuple val(meta), path("${meta.tag}.gff3.tar.gz"),    emit: tar_gff3_annotations, optional: true // bundled
    // tuple val(meta), path("${meta.tag}.json.tar.gz"),    emit: tar_json_annotations, optional: true // bundled
    tuple val(meta), path("*.gff3"),                     emit: gff3_annotations,     optional: true // unbundled
    tuple val(meta), path("SORFs_bakta/*.json"),         emit: json_annotations,     optional: true // unbundled
    tuple val(meta), path("${meta.tag}.sorf-extra.pkl"), emit: pkl_annotations
    tuple val(meta), path("${meta.tag}.json.tar"), path("${meta.tag}.genomes.tar"),      emit: tar_archives_for_bakpack, optional: true // bundled -> bakpack


    script:
    def compliant = params.compliant ? "--compliant" : ""
    def individual_pickles = meta.asm_ids.collect { asm_id -> "SORFs_bakta/${asm_id}.sorf-extra.pkl" }.join(' ')
    def individual_jsons = meta.asm_ids.collect { asm_id -> "SORFs_bakta/${asm_id}.json" }.join(' ')
    def individual_gff3s = meta.asm_ids.collect { asm_id -> "${asm_id}.gff3" }.join(' ')
    def individual_fastas = meta.asm_ids.collect { asm_id -> "${asm_id}.fa" }.join(' ')
    def rm_individual_gff3s = params.save_intermediate ? "" : "&& rm -rf ${individual_gff3s}"
    def rm_individual_jsons = params.save_intermediate ? "" : "&& rm -rf ${individual_jsons}"
    def rm_individual_pickles = params.save_intermediate ? "" : "&& rm -rf ${individual_pickles}"
    def rm_additional_intermediates = params.save_intermediate ? "" : "&& rm -rf cds_unbatched/*.cds.pkl rna_unbatched/*.rna.pkl SORFs_bakta/*.log"
    """
    manage_pkls.py unbatch \\
        --input ${cds_pkl} \\
        --out-dir cds_unbatched \\
        --suffix .cds.pkl

    manage_pkls.py unbatch \\
        --input ${rna_pkl} \\
        --out-dir rna_unbatched \\
        --suffix .rna.pkl

    # Loop through assemblies in batch running bakta on each
    for asm in ${assemblies}; do
        id=\$(basename "\$asm" | sed -E 's/(\\.[^.]+)+\$//')
        bakta --db ${bakta_db} --sorf-extra  \\
            --cds-pickle cds_unbatched/\${id}.cds.pkl \\
            --rna-pickle rna_unbatched/\${id}.rna.pkl \\
            --output SORFs_bakta/  \\
            --prefix \${id} \\
            --threads ${task.cpus} \\
            --force \\
            ${compliant} \\
            \${asm}
        mv SORFs_bakta/\${id}.sorf-extra.gff3 \${id}.gff3
    done

    # Concatenate to batch-level files
    manage_pkls.py batch \\
        --sample-ids ${meta.asm_ids.join(',')} \\
        --output ${meta.tag}.sorf-extra.pkl  \\
        ${individual_pickles}

    # Bundle outputs and remove intermediates
    if ( "${params.bundle_gff3}" == "true" ) ; then
        tar -czf ${meta.tag}.gff3.tar.gz ${individual_gff3s}  && \
        tar -cf ${meta.tag}.json.tar -C SORFs_bakta ${individual_jsons} && \
        tar -chf ${meta.tag}.genomes.tar ${individual_fastas} \
            ${rm_individual_gff3s} \
            ${rm_individual_jsons} \
            ${rm_individual_pickles} \
            ${rm_additional_intermediates}
    fi
    """
}