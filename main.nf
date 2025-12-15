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
    cds_dir = Channel.of("/path/to/sample1.fa.cds-only.faa", "/path/to/sample1.fa.cds-only.faa") // provide individual FASTA files with found CDS

    CLUSTER_PROTEOME(cds_dir)
    CLUSTER_PROTEOME
        .out
        .map { all_seqs, clustering_tsv, rep_seq -> rep_seq }
        .set { clustered_proteins_ch }
}