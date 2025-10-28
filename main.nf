#!/usr/bin/env nextflow
// Copyright (C) 2024 Genome Research Ltd.

/*
========================================================================================
    HELP
========================================================================================
*/

def logo = NextflowTool.logo(workflow, params.monochrome_logs)

log.info logo

NextflowTool.commandLineParams(workflow.commandLine, log, params.monochrome_logs)


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

// include { FIND_CDSS } from './subworkflows/find_cdss.nf'
include { ANNOTATE_PROTEINS } from './subworkflows/annotate_proteins.nf'
// include { CLUSTER_PROTEOME } from './subworkflows/proteome_clustering.nf'
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
    
    infiles = Channel.fromPath("anonymised_proteins_for_debug/some_proteins.faa")
        // .view() // DEBUG

    // DEBUG test annotation process
    ANNOTATE_PROTEINS(infiles)

    // GENERATE_PANGENOME(infiles)

    // GENERATE_PANGENOME.out
    //     .set { pangenome_index }

    // ANNOTATE_USING_PANGENOME(pangenome_index)
    
    
    // FIND_CDSS(infiles) | CLUSTER_PROTEOME | GENERATE_PANGENOME | ANNOTATE_USING_PANGENOME
}