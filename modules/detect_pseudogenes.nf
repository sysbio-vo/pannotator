process DETECT_PSEUDOGENES {
    tag { meta.tag }
    label "cds_pseudogenes"
    label 'bakta'
    publishDir "${params.outdir}/predicted_cds_with_pseudogenes", enabled: params.save_intermediate, mode: 'copy'

    input:
    tuple val(meta), path(batch_pickle_manifest)
    path(bakta_db)

    output:
    tuple val(meta), path("cds_with_pseudogenes/*with_pseudogenes.pkl")

    script:
    // output_prefix = assembly.getBaseName()
    // bakta_db_arg = params.bakta_db ? "--db ${params.bakta_db}" : ""
    out_pickle_dir = "cds_with_pseudogenes"
    """
    # Loop through assemblies in batch running bakta on each
    bakta_pseudo_bulk \
      --db ${bakta_db} \
      --output ${out_pickle_dir} \
      --prefix pseudogenes \
      --batch_pickle_manifest ${batch_pickle_manifest}
    """
}