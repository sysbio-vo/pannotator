#!/usr/bin/env nextflow
// Copyright (C) 2024 Genome Research Ltd.

/*
========================================================================================
    HELP
========================================================================================
*/

// def logo = nextflowtool.NextflowTool.logo(workflow, params.monochrome_logs)

// log.info logo

// nextflowtool.NextflowTool.commandLineParams(workflow.commandLine, log, params.monochrome_logs)


def printHelp() {
    nextflowtool.NextflowTool.help_message("${workflow.ProjectDir}/schema.json", 
                               [],
    params.monochrome_logs, log)
}

def sampleIdFromName = {name -> name.replaceFirst(~/(\.[^\.]+)+$/, '')}


// if a Bakta DB doesn't exist, this value will be null
def parsed_bakta_db_type = utils.PannotatorUtils.get_bakta_db_type("${params.bakta_db}")
// if the DB doesn't exist, it will be downloaded using provided type
def bakta_db_type = parsed_bakta_db_type ? parsed_bakta_db_type : params.bakta_db_type 

/*
========================================================================================
    IMPORT MODULES/SUBWORKFLOWS
========================================================================================
*/

include { BATCHER as FIND_CDSS_BATCHER;
          BATCHER as FIND_RNAS_BATCHER;
          BATCHER as SORF_EXTRA_BATCHER } from './subworkflows/helpers.nf'
include { FIND_CDSS } from './subworkflows/find_cdss.nf'
include { ANNOTATE_PROTEINS; ANNOTATE_WITH_AUX_DB } from './subworkflows/annotate_proteins.nf'
include { CLUSTER_PROTEOME } from './subworkflows/proteome_clustering.nf'
include { MERGE_ANNOTATIONS } from './modules/merge_annotations.nf'
include { DETECT_PSEUDOGENES_OPTIONAL as DETECT_PSEUDOGENES } from './subworkflows/detect_pseudogenes.nf'
include { FIND_RNAS } from './modules/find_rnas.nf'
include { DOWNLOAD_BAKTA_DB } from './modules/helpers.nf'
include { SORF_EXTRA } from './modules/find_sorf_extra.nf'
include { EXTEND_OR_GENERATE_AUXILIARY_DB } from './modules/generate_auxiliary_db.nf'
include { EXTEND_ANNOTATIONS } from './modules/extend_annotations.nf'


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

    infiles = Channel.fromPath("${params.indir}/*${params.infile_extension}")

    // Load assemblies into batches and track this metadata
    find_cds_batches_and_bakta_db = FIND_CDSS_BATCHER(infiles, bakta_db, params.find_cds_buffer_size, "find_cds")
    find_rna_batches_and_bakta_db = FIND_RNAS_BATCHER(infiles, bakta_db, params.find_rnas_buffer_size, "find_rnas")
    sorf_extra_batches_and_bakta_db = SORF_EXTRA_BATCHER(infiles, bakta_db, params.sorf_extra_buffer_size, "sorf_extra")

    //ch_asm_by_id = infiles.map { asm -> tuple(sampleIdFromName(asm.name), asm) }

    //-----------------------------
    // CDS prediction
    //-----------------------------
    cds_outputs = FIND_CDSS(batches_and_bakta_db)
    batch_cds_faas = cds_outputs.map { meta, faa, pkl -> faa } // no meta as all inputs pool in CLUSTER_PROTEOME
    batch_cds_pkls = cds_outputs.map { meta, faa, pkl -> tuple(meta, pkl) }

    //-----------------------------
    // Cluster + annotate
    //-----------------------------
    CLUSTER_PROTEOME(batch_cds_faas)
    CLUSTER_PROTEOME.out.set { cluster_out_ch }  

    // TODO: use multiMap (https://docs.seqera.io/nextflow/reference/operator#multimap)
    rep_proteins_ch = cluster_out_ch.map { all_seqs, clustering_tsv, rep_seq -> rep_seq }
    clustering_tsv_ch = cluster_out_ch.map { all_seqs, clustering_tsv, rep_seq -> clustering_tsv }
    all_seqs_ch = cluster_out_ch.map { all_seqs, clustering_tsv, rep_seq -> all_seqs }   

    rep_proteins_ch
        .combine(bakta_db)
        .set { rep_proteins_and_dbs }

    hmm_ch = params.user_hmms ? Channel.of(file(params.user_hmms)) : []
    prots_ch = params.user_proteins ? Channel.of(file(params.user_proteins)) : []

    // if an auxiliary db path is provided, utilise it
    // to propagate existing annotations into the annotation JSON
    if ( params.auxiliary_db && file(params.auxiliary_db).exists() ) {

        auxiliary_db = Channel.of(file(params.auxiliary_db))

        rep_proteins_and_dbs
            .combine(auxiliary_db)
            .set { rep_proteins_and_dbs }

        ANNOTATE_WITH_AUX_DB(rep_proteins_and_dbs, hmm_ch, prots_ch)
        ANNOTATE_WITH_AUX_DB.out.bulk_annotations
            .set { bulk_annotations }
    } else {
        ANNOTATE_PROTEINS(rep_proteins_and_dbs, hmm_ch, prots_ch)
        ANNOTATE_PROTEINS.out.bulk_annotations
            .set { bulk_annotations }
    }

    //-----------------------------
    // Extend annotations to cluster members
    // (for non-identical clustering)
    //-----------------------------

    // TODO: this is unreliable
    // consider switching to explicit parameter definition in the config
    if( (params.mmseqs_args ?: '') != '--min-seq-id 1.0 -c 1.0 --alignment-mode 3' ) {
        EXTEND_ANNOTATIONS(
            clustering_tsv_ch,
            all_seqs_ch,
            bulk_annotations
        )
        bulk_ann_final_ch = EXTEND_ANNOTATIONS.out.bulk_annotations_extended
    } else {
        bulk_ann_final_ch = bulk_annotations
    }

    //-----------------------------
    // Merge annotations
    //-----------------------------
    
    // NOTE: cache is not utilised if channel values are collected in a different order
    // TODO: sort collected values in cds_pkl_list_ch?
    MERGE_ANNOTATIONS(
        batch_cds_pkls,
        bulk_ann_final_ch
    )

    DETECT_PSEUDOGENES(MERGE_ANNOTATIONS.out.annotated_pickles, bakta_db, bakta_db_type)
    DETECT_PSEUDOGENES.out.annotated_samples_updated
        .set { ch_cds_annot_pkl }

    //-----------------------------
    // RNA prediction
    //-----------------------------
    batch_rna_pkls = FIND_RNAS(batches_and_bakta_db)
    
    //-----------------------------
    // SORF extra search
    //-----------------------------
    ch_cds_keyed = ch_cds_annot_pkl.map { p -> tuple(sampleIdFromName(p.name), p) }
    ch_rna_keyed = ch_rna_pkl.map { p -> tuple(sampleIdFromName(p.name), p) }

    ch_sorf_in = ch_cds_keyed
        .join(ch_rna_keyed)
        .join(ch_asm)
        .map { sid, cds_pkl, rna_pkl, asm -> tuple(sid, asm, cds_pkl, rna_pkl) }
        .combine(bakta_db)

    SORF_EXTRA(ch_sorf_in)

    // TODO: refactor branching
    if ( params.auxiliary_db && (!file(params.auxiliary_db).exists() || params.extend_auxdb) ) {
        SORF_EXTRA
            .out
            .gff3_annotations
            .map { sample_id, anno_gff3, anno_pkl -> anno_pkl }
            .collect()
            .set { final_pkl_anno }

        auxiliary_db = Channel.of(file(params.auxiliary_db)).combine(bulk_annotations)

        EXTEND_OR_GENERATE_AUXILIARY_DB(final_pkl_anno, ch_cds_annot_pkl, auxiliary_db)
    }
}
