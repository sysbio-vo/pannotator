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
include { BUILD_COORDS_INDEX_WF } from './subworkflows/build_coords_index_wf.nf'

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
    infiles = Channel.fromPath("${params.indir}/*")
        // .take( 10 ) // DEBUG
        // .view() // DEBUG
    outdir = file(params.outdir)

    
    cdss_dir = FIND_CDSS(infiles, Channel.value(outdir))
    
    cds_dir = FIND_CDSS.out
                .collect()
                .map { outdir.resolve('CDSS_bakta') }
    BUILD_COORDS_INDEX_WF(cds_dir, Channel.value(outdir))

    // GENERATE_PANGENOME(infiles)

    // GENERATE_PANGENOME.out
    //     .set { pangenome_index }

    // ANNOTATE_USING_PANGENOME(pangenome_index)
    
    
    // FIND_CDSS(infiles) | CLUSTER_PROTEOME | GENERATE_PANGENOME | ANNOTATE_USING_PANGENOME
}