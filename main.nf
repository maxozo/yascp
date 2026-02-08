#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/
include { YASCP } from "$projectDir/workflows/yascp"
include { RETRIEVE_RECOURSES; RETRIEVE_RECOURSES_TEST_DATASET } from "$projectDir/modules/local/retrieve_resources/retrieve_resources"
include { RSYNC_RESULTS_REMOVE_WORK_DIR} from "$projectDir/modules/local/rsync_results_remove_work_dir/main"
include { CELLTYPE} from "$projectDir/subworkflows/celltype"
include { QC_AND_INTEGRATION } from "$projectDir/subworkflows/qc_and_integration"
include { DUMMY_FILTERED_CHANNEL } from "$projectDir/modules/local/merge_samples/functions"
include { MATCH_GT_VIREO } from "$projectDir/modules/local/genotypes/main"
include { MULTIPLET } from "$projectDir/subworkflows/doublet_detection"
include { ENHANCE_STATS_GT_MATCH } from "$projectDir/modules/local/genotypes/main"

/**
 * Helper to parse VCF TSVs consistently
 */
def parseVcfTsv(tsv_path) {
    return Channel.fromPath(tsv_path, followLinks: true, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row -> 
            tuple(row.label, file(row.vcf_file_path), file("${row.vcf_file_path}.csi")) 
        }
}

/*
========================================================================================
    MAIN ENTRY POINT
========================================================================================
*/

workflow MAIN {
    // 1. Set up output directory
    def outdir_path = params.outdir ?: "${launchDir}/${outdir}"
    out_ch = Channel.fromPath(outdir_path, checkIfExists: true)

    // 2. Define Input Channels
    if (params.profile == 'test_full') {
        RETRIEVE_RECOURSES_TEST_DATASET(out_ch)
        input_channel = RETRIEVE_RECOURSES_TEST_DATASET.out.input_channel
        vcf_inputs = RETRIEVE_RECOURSES_TEST_DATASET.out.vcf_inputs
            .splitCsv(header: true, sep: '\t')
            .map { row -> tuple(row.label, file(row.vcf_file_path), file("${row.vcf_file_path}.csi")) }
    } else {
        input_channel = Channel.fromPath(params.input_data_table, followLinks: true, checkIfExists: true)
        
        vcf_inputs = params.genotype_input.run_with_genotype_input
            ? parseVcfTsv(params.genotype_input.tsv_donor_panel_vcfs)
            : Channel.empty()
    }

    // 3. Track inputs and trigger main YASCP workflow
    input_channel.collectFile(name: "${params.outdir}/yascp_inputs.tsv")
    
    // Call the workflow included from workflows/yascp.nf
    YASCP('default', input_channel, vcf_inputs)
}

workflow {
    MAIN()
}

/*
========================================================================================
    NAMED ENTRY POINTS (Reruns/Specific Tasks)
========================================================================================
*/

workflow JUST_CELLTYPES {
    CELLTYPE(Channel.fromPath(params.file__anndata_merged), 'celltype_mode')
}

workflow JUST_DOUBLETS {
    MULTIPLET(Channel.fromPath(params.file__anndata_merged), 'doublet_mode')
}

workflow GT_MATCH {
    parseVcfTsv(params.genotype_input.tsv_donor_panel_vcfs)
        .map { label, vcf, csi -> 
            tuple(params.file_name, file(params.vcf), file("${params.vcf}.tbi"), label, vcf, csi) 
        }
        .set { gt_match_input }

    MATCH_GT_VIREO(gt_match_input)
    ENHANCE_STATS_GT_MATCH(MATCH_GT_VIREO.out.donor_match_table_with_pool_id, params.input_data_table)
}

workflow WORK_DIR_REMOVAL {
    def outdir_path = params.outdir ?: "${launchDir}/${outdir}"
    RSYNC_RESULTS_REMOVE_WORK_DIR(Channel.fromPath(outdir_path, checkIfExists: true), params.tmpdir)
}

/*
========================================================================================
    COMPLETION HANDLER
========================================================================================
*/

workflow.onComplete {
    log.info """
    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    """
    
    if (workflow.success && params.remove_work_dir) {
        log.info "Cleaning up work directory: ${params.tmpdir}"
        "bash ${projectDir}/bin/del_work_dirs.sh ${params.tmpdir}".execute().waitFor()
    }
}