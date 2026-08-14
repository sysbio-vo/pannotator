process SORF_EXTRA {
    tag { sample_id }
    label "sorf_extra_search"
    label 'bakta'

    // TODO: allow specifying output path for the gff3 file in Bakta to avoid doing this
    publishDir (
        "${params.outdir}/final_annotations", 
        pattern: "*.sorf-extra.gff3",
        mode: 'copy', 
        saveAs: { filename -> 
            def name = file(filename).name
            if (name.endsWith(".sorf-extra.gff3")) {
                def base = name.replaceFirst(/\.sorf-extra\.gff3$/, '')
                return "${base}.gff3"
            }
            return null
        }, enabled: { ${params.bundle_gff3} != true }
    )
    publishDir (
        "${params.outdir}/final_annotations", 
        pattern: "*.gff3.tar.gz",
        mode: 'copy', 
        enabled: params.bundle_gff3
    )
    publishDir (
        "${params.outdir}/final_annotations", 
        pattern: "*.sorf-extra.pkl",
        mode: 'copy', 
        enabled: params.save_intermediate
    )

    input:
    tuple val(meta), path(assemblies), path(cds_pkl), path(rna_pkl), path(bakta_db)

    output:
    tuple val(meta), path("${meta.tag}.gff3.tar.gz"), emit: gff3_annotations
    tuple val(meta), path("${meta.tag}.sorf-extra.pkl"), emit: pkl_annotations


    script:
    def compliant = params.compliant ? "--compliant" : ""
    def individual_pickles = meta.asm_ids.collect { asm_id -> "SORFS_bakta/${asm_id}.sorf-extra.pkl" }.join(' ')
    def individual_gff3s = meta.asm_ids.collect { asm_id -> "SORFS_bakta/${asm_id}.sorf-extra.gff3" }.join(' ')
    """
    # Loop through assemblies in batch running bakta on each
    for asm in ${assemblies}; do
        id=\$(basename "\$asm" | sed -E 's/(\\.[^.]+)+\$//')
        bakta --db ${bakta_db} --sorf-extra  \
            --cds-pickle ${cds_pkl} \
            --rna-pickle ${rna_pkl} \
            --output SORFS_bakta/  \
            --prefix ${id} \
            --threads ${task.cpus} \
            ${compliant} \
            ${asm}
    done

    # Concatenate to batch-level files
    merge_pickles.py \\
        --assembly_ids ${meta.asm_ids.join(',')} \\
        --out ${meta.tag}.sorf-extra.pkl \\
        ${individual_pickles} \
        && rm -rf ${individual_pickles}

    if ( "${params.bundle_gff3}" == "true" ) ; then
        tar -czf ${meta.tag}.gff3.tar.gz ${individual_gff3s} \
            --transform "s|SORFS_bakta/||" \
            && rm -rf ${individual_gff3s}
    fi
    """
}