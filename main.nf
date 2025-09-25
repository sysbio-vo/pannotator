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

include { BUILD_PANGENOME } from './modules/generate_pangenome.nf'

include { FIND_CDSS } from './subworkflows/find_cds.nf'
include { CLUSTER_PROTEOME } from './subworkflows/proteome_clustering.nf'
include { ANNOTATE_USING_PANGENOME } from './subworkflows/pangenome_annotation.nf'



/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

workflow PANNOTATE {
    if (params.help) {
        printHelp()
        exit 0
    }
    
    infiles = Channel.fromPath('${params.indir}/*')

    GENERATE_PANGENOME(infiles)

    GENERATE_PANGENOME.out
        .set { pangenome_index }

    ANNOTATE_USING_PANGENOME(pangenome_index)
    
    
    // FIND_CDSS(infiles) | CLUSTER_PROTEOME | GENERATE_PANGENOME | ANNOTATE_USING_PANGENOME
}