process SORF_EXTRA {
    tag { sample_id }
    label "sorf_extra_search"
    label 'bakta'

    // TODO: allow specifying output path for the gff3 file in Bakta to avoid doing this
    publishDir (
        params.outdir, 
        mode: 'copy', 
        saveAs: { filename -> 
            def name = file(filename).name
            if (name.endsWith(".sorf-extra.gff3")) {
                def base = name.replaceFirst(/\.sorf-extra\.gff3$/, '')
                return "${base}.gff3"
            }
            // // TODO: save final pickle objects 
            // else if (name.endsWith(".sorf-extra.pkl")) {
            //     return name
            // }
            return null
        }
    )

    input:
    tuple val(sample_id), path(assembly), path(cds_pkl), path(rna_pkl), path(bakta_db)

    output:
    tuple val(sample_id), path("${sample_id}/${sample_id}.sorf-extra.gff3"), path("${sample_id}/${sample_id}.sorf-extra.pkl"), emit: gff3_annotations

    script:
    def compliant = params.compliant ? "--compliant" : ""
    """
    bakta --db ${bakta_db} --sorf-extra  \
        --cds-pickle ${cds_pkl} \
        --rna-pickle ${rna_pkl} \
        --output ${sample_id}  \
        --prefix ${sample_id} \
        --threads ${task.cpus} \
        ${compliant} \
        ${assembly}
    """
}