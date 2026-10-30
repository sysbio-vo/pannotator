process SORF_EXTRA {
    tag { meta.tag }
    label "sorf_extra_search"
    label 'bakta'

    publishDir (
        "${params.outdir}/final_annotations",
        pattern: "*.gff3.tar.gz",
        enabled: params.bundle_type in ["gff3", "both"],
        mode: 'copy'
    )
    publishDir (
        "${params.outdir}/final_annotations",
        pattern: "*.gff3",
        enabled: params.bundle_type in ["null", null],
        mode: 'copy'
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
    tuple val(meta), path("${meta.tag}.gff3.tar.gz"),                                       emit: tar_gff3_annotations, optional: true // bundled gff3
    tuple val(meta), path("${meta.tag}.json.tar.xz"), path("${meta.tag}.genomes.tar.xz"),   emit: tar_json_annotations, optional:true // bundled json
    tuple val(meta), path("*.gff3"),                                                        emit: gff3_annotations,     optional: true // unbundled
    tuple val(meta), path("${meta.tag}.sorf-extra.pkl"),                                    emit: pkl_annotations


    script:
    def compliant = params.compliant ? "--compliant" : ""

    // bundle logic, can be bakpack, gff3, both, or null (individual gff3 files by default)
    // If bakpack is selected bundle together json files to be passed to the bakpack process
    def bundle_type = params.bundle_type?.toString()?.toLowerCase()
    if (!(bundle_type in ['bakpack', 'gff3', 'both', null, 'null'])) {
        error "Unknown --bundle_type '${params.bundle_type}' (expected bakpack, gff3, both or null)"
    }
    def tar_gff3 = bundle_type in ['gff3', 'both']
    def tar_json = bundle_type in ['bakpack', 'both']
    def tar_fa = bundle_type in ['bakpack', 'both']

    // Collect required bundles
    def individual_jsons = meta.asm_ids.collect { asm_id -> "${asm_id}.json" }.join(' ')
    def individual_gff3s = meta.asm_ids.collect { asm_id -> "SORFs_bakta/${asm_id}.sorf-extra.gff3" }.join(' ')
    def individual_fas = meta.asm_ids.collect { asm_id -> "${asm_id}.fa" }.join(' ')
    def individual_pickles = meta.asm_ids.collect { asm_id -> "SORFs_bakta/${asm_id}.sorf-extra.pkl" }.join(' ')

    // process file clean up commands
    def rm_individual_gff3s = params.save_intermediate ? "" : "&& rm -rf ${individual_gff3s}"
    def rm_individual_jsons = params.save_intermediate ? "" : "&& rm -rf ${individual_jsons}"
    def rm_individual_fas = params.save_intermediate ? "" : "&& rm -rf ${individual_fas}"

    // perform bundles
    def gff3_cmd = tar_gff3 ? "tar -czf ${meta.tag}.gff3.tar.gz ${individual_gff3s} ${rm_individual_gff3s}" : ""
    def json_cmd = tar_json ? "tar -C SORFs_bakta -cf -  ${individual_jsons} | xz -T ${task.cpus} -9e > ${meta.tag}.json.tar.xz ${rm_individual_jsons}" : ""
    def fna_cmd = tar_fa ? "tar -chf - ${individual_fas} | xz -T ${task.cpus} -9e > ${meta.tag}.genomes.tar.xz ${rm_individual_fas}" : ""

    // tar -cf - *.fna | xz -T 1 -9e -c > genomes.tar.xz
    // tar -cf - *.json | xz -T 1 -9e -c > annotations.tar.xz
    // bakpack build   --annotations annotations.tar.xz   --genomes genomes.tar.xz   --xz-threads 1   --output annotations.bakpack
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

    # Bundle annotations
    ${gff3_cmd}
    
    export LC_ALL=C
    ${json_cmd}
    ${fna_cmd}
    """
}