#!/usr/bin/env nextflow
// Copyright (C) 2024 Genome Research Ltd.

/*
========================================================================================
    HELP
========================================================================================
*/

// def logo = NextflowTool.logo(workflow, params.monochrome_logs)

// log.info logo

// NextflowTool.commandLineParams(workflow.commandLine, log, params.monochrome_logs)


def printHelp() {
    NextflowTool.help_message("${workflow.ProjectDir}/schema.json", 
                               [],
    params.monochrome_logs, log)
}

/*
========================================================================================
    IMPORT MODULES/SUBWORKFLOWS
========================================================================================
*/

// include { BUILD_PANGENOME } from './modules/generate_pangenome.nf'

include { FIND_CDSS } from './subworkflows/find_cdss.nf'
include { ANNOTATE_PROTEINS } from './subworkflows/annotate_proteins.nf'
include { BUILD_COORDS_INDEX_WF } from './subworkflows/build_coords_index_wf.nf'
include { CLUSTER_PROTEOME } from './subworkflows/proteome_clustering.nf'
include { MERGE_ANNOTATIONS } from './modules/merge_annotations.nf'
include { DETECT_PSEUDOGENES } from './subworkflows/detect_pseudogenes.nf'
include { FIND_RNAS } from './modules/find_rnas.nf'
include { DOWNLOAD_BAKTA_DB } from './modules/helpers.nf'

// include { ANNOTATE_USING_PANGENOME } from './subworkflows/pangenome_annotation.nf'



/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow {
    if (params.help) {
        printHelp()
        exit 0
    }

    // TODO: even wich `cache` directive set to false
    // the database is stored twice - in the workdir and publishdir
    // use the bakta_db config parameter as input for subsequent processes
    // instead of the output of the DOWNLOAD_BAKTA_DB process
    // and set publishDir move to `move` in that process
    if ( file(params.bakta_db).exists() ) {
        bakta_db = Channel.of(file(params.bakta_db))
    } else {
        println "Downloading bakta db to ${params.bakta_db}"
        bakta_db = DOWNLOAD_BAKTA_DB(params.bakta_db_type)
    }

    infiles = Channel.fromPath("${params.indir}/*")
        .take( 3 ) // DEBUG
        // .view() // DEBUG

    infiles
        .combine(bakta_db)
        .set { infiles_and_bakta_db }

    cds_outputs = FIND_CDSS(infiles_and_bakta_db)

    all_cds_outputs = cds_outputs.collect()

    cds_pkl_files = all_cds_outputs
        .flatten()
        .filter { it.name.endsWith('.pkl') }
        .collect()

    // BUILD_COORDS_INDEX_WF(cds_pkl_files) 
    // cds_outputs.view() // DEBUG

    CLUSTER_PROTEOME(cds_outputs)
    CLUSTER_PROTEOME
        .out
        .map { all_seqs, clustering_tsv, rep_seq -> rep_seq }
        .set { rep_proteins_ch }

    rep_proteins_ch
        .combine(bakta_db)
        .set { rep_proteins_and_bakta_db }

    ANNOTATE_PROTEINS(rep_proteins_and_bakta_db)

    rna_outputs = FIND_RNAS(infiles)
    all_rna_outputs = rna_outputs.collect()
    rna_pickle_files = all_rna_outputs
        .flatten()
        .collect()

    MERGE_ANNOTATIONS(
        cds_pkl_files,
        ANNOTATE_PROTEINS.out.bulk_annotations
    )

    // predict pseudogenes using annotated pickle objects
    MERGE_ANNOTATIONS.out
        .flatten()
        .map { it -> "${it}" } // TODO: is this crutch REALLY neccessary to collect paths to files in a txt files instead of their contents? 
        .collectFile( name: 'annotated_cds_manifest.txt', newLine: true )
        .set { manifest_file }
    
    manifest_file
        .combine(bakta_db)
        .set { manifest_file_and_bakta_db }


    DETECT_PSEUDOGENES(manifest_file_and_bakta_db)
}