process SORF_EXTRA {
    tag { meta.tag }
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
        }, enabled: { params.bundle_gff3 != true }
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
    tuple val(meta), path("${meta.tag}.gff3.tar.gz"), emit: tar_gff3_annotations, optional: true
    tuple val(meta), path("${meta.tag}.gff3"), emit: gff3_annotations, optional: true
    tuple val(meta), path("${meta.tag}.sorf-extra.pkl"), emit: pkl_annotations


    script:
    def compliant = params.compliant ? "--compliant" : ""
    def individual_pickles = meta.asm_ids.collect { asm_id -> "SORFs_bakta/${asm_id}.sorf-extra.pkl" }.join(' ')
    def individual_gff3s = meta.asm_ids.collect { asm_id -> "SORFs_bakta/${asm_id}.sorf-extra.gff3" }.join(' ')
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
    done

    # Concatenate to batch-level files
    manage_pkls.py batch \\
        --sample-ids ${meta.asm_ids.join(',')} \\
        --output ${meta.tag}.sorf-extra.pkl  \\
        ${individual_pickles}

    if ( "${params.bundle_gff3}" == "true" ) ; then
        tar -czf ${meta.tag}.gff3.tar.gz ${individual_gff3s} \\
            --transform "s|SORFs_bakta/||" \\
            && rm -rf ${individual_gff3s}
    fi
    """
}