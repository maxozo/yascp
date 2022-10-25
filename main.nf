#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/yascp
========================================================================================
    Github : https://github.com/nf-core/yascp
    Website: https://nf-co.re/yascp
    Slack  : https://nfcore.slack.com/channels/yascp
----------------------------------------------------------------------------------------
*/


nextflow.enable.dsl = 2

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/


/*
========================================================================================
  NAMED WORKFLOWS FOR TESTS
========================================================================================
*/

//include { TEST_MATCH_GENOTYPES } from './tests/test_genotypes'
//include { TEST_GT_ASSIGN } from './tests/test_gt_assign'
//include { TEST_GTCHECK } from './tests/test_gtcheck'
//include { TEST_MATCH_GT_VIREO } from './tests/test_genotypes'
//include { TEST_SPLIT_BAM_PER_DONOR } from './tests/test_bam_per_donor'
//include { TEST_ENCRYPT_DIR} from './tests/test_encryption'
include { TEST_SUBSET_GENOTYPES } from "$projectDir/tests/test_subset_genotypes"
include { match_genotypes } from "$projectDir/subworkflows/match_genotypes"
include { metadata_posthoc;replace_donors_posthoc } from "$projectDir/modules/local/report_update/main"
include {data_handover} from "$projectDir/subworkflows/data_handover"
include { CREATE_ARTIFICIAL_BAM_CHANNEL } from "$projectDir/modules/local/create_artificial_bam_channel/main"
include { TRANSFER;SUMMARY_STATISTICS_PLOTS } from "$projectDir/modules/nf-core/modules/summary_statistics_plots/main"
include {SUBSET_WORKF} from "$projectDir/modules/nf-core/modules/subset_genotype/main"
include {VIREO} from "$projectDir/modules/nf-core/modules/vireo/main"
include {capture_cellbender_files} from "$projectDir/modules/nf-core/modules/cellbender/functions"
include { DECONV_INPUTS } from "$projectDir/subworkflows/prepare_inputs/deconvolution_inputs"
include { prepare_inputs } from "$projectDir/subworkflows/prepare_inputs"
include {MULTIPLET} from "$projectDir/modules/nf-core/modules/multiplet/main"
include { SPLIT_DONOR_H5AD } from "$projectDir/modules/nf-core/modules/split_donor_h5ad/main"
include {REPLACE_GT_DONOR_ID2 } from "$projectDir/modules/nf-core/modules/genotypes/main"

workflow NF_CORE_TEST {
  //println "**** running NF_CORE_TEST::TEST_MATCH_GENOTYPES"
  //TEST_MATCH_GENOTYPES()
  //println "**** running NF_CORE_TEST::TEST_GT_ASSIGN"
  //TEST_GT_ASSIGN()
  //println "**** running NF_CORE_TEST::TEST_GTCHECK"
  //TEST_GTCHECK()
  //println "**** running NF_CORE_TEST::TEST_MATCH_GT_VIREO"
  //TEST_MATCH_GT_VIREO()
  //println "**** running NF_CORE_TEST::TEST_SPLIT_BAM_PER_DONOR"
  //TEST_SPLIT_BAM_PER_DONOR()
  //println "**** running NF_CORE_TEST::TEST_ENCRYPT_DIR"
  //TEST_ENCRYPT_DIR()
    println "**** running NF_CORE_TEST::TEST_SUBSET_GENOTYPES"
    //   TEST_SUBSET_GENOTYPES()
    input_channel = Channel.fromPath(params.input_data_table, followLinks: true, checkIfExists: true)
    
    prepare_inputs(input_channel)

    capture_cellbender_files(params.cellbender_location,"${params.output_dir}/nf-preprocessing")
    DECONV_INPUTS(capture_cellbender_files.out.celbender_path,prepare_inputs)
    channel__file_paths_10x = DECONV_INPUTS.out.channel__file_paths_10x
    ch_experiment_filth5= DECONV_INPUTS.out.ch_experiment_filth5
    MULTIPLET(
        params.output_dir,
        channel__file_paths_10x,
        params.sample_qc.cell_filters.filter_multiplets.expected_multiplet_rate,
        params.sample_qc.cell_filters.filter_multiplets.n_simulated_multiplet,
        params.sample_qc.cell_filters.filter_multiplets.multiplet_threshold_method,
        params.sample_qc.cell_filters.filter_multiplets.scale_log10
    )

    log.info "#### running DECONVOLUTION workflow #####"
    if (params.run_with_genotype_input) {
        // if (params.genotype_input.subset_genotypes){
            log.info "---We are subsetting genotypes----"
            // We have to produce a single vcf file for each individual pool.
            // Therefore we create 2 channels:
            // 1) All the expected vcf ids listed in the donor table
            Channel.fromPath(params.input_data_table,      
            followLinks: true,
            checkIfExists: true
            ).splitCsv(header: true, sep: '\t').map { row -> tuple(row.experiment_id, row.donor_vcf_ids) }
            .set { donors_in_pools }
            
            // 2) All the vcfs provided to us. 
            Channel.fromPath(
            params.genotype_input.tsv_donor_panel_vcfs,
            followLinks: true,
            checkIfExists: true
            ).splitCsv(header: true, sep: '\t')
            .map { row -> tuple(row.label, file(row.vcf_file_path), file("${row.vcf_file_path}.csi")) }
            .set { ch_ref_vcf }

            // This will subsequently result in a joint vcf file for all the cohorts listed for each of the pools that can be used in VIREO and/or GT matching algorythm.
            SUBSET_WORKF(ch_ref_vcf,donors_in_pools)
            merged_expected_genotypes = SUBSET_WORKF.out.merged_expected_genotypes
            // merged_expected_genotypes.view()
            // SUBSET_GENOTYPE(ch_experiment_donorsvcf_donorslist.map { experiment, donorsvcf,donortbi, donorslist -> tuple(experiment,
            //                 file(donorsvcf),file(donortbi),
            //                 donorslist)})
        // }
    }

    channel_input_data_table = Channel.fromPath(params.input_data_table, followLinks: true, checkIfExists: true)
        
    if (params.existing_cellsnp != ''){
        expl1 = Channel.fromPath( "${params.existing_cellsnp}/*/*.vcf.gz")
        expl1.map{row -> tuple("${row[-2]}".replaceAll('cellsnp_',''), "${row}".replaceAll('/cellSNP.base.vcf.gz',''))}.set{cellsnp_output_dir}
        
    }else{
        log.inf('Running CELLSNP')
        channel_input_data_table
            .splitCsv(header: true, sep: params.input_tables_column_delimiter)
            .map{row->tuple(row.experiment_id, "${row.data_path_10x_format}/possorted_genome_bam.bam" ,row.data_path_10x_format+'/filtered_feature_bc_matrix/barcodes.tsv.gz')}
            .set{pre_ch_experiment_bam_barcodes}
        pre_ch_experiment_bam_barcodes
            .map { a,b,c -> tuple(a, file(b), file("${b}.bai"), file(c))}
            .set {ch_experiment_bam_bai_barcodes}

        CELLSNP(ch_experiment_bam_bai_barcodes,
            Channel.fromPath(params.cellsnp.vcf_candidate_snps).collect())
        cellsnp_output_dir = CELLSNP.out.cellsnp_output_dir
    }

    channel_input_data_table
        .splitCsv(header: true, sep: params.input_tables_column_delimiter)
        .map{row->tuple(row.experiment_id, row.n_pooled)}
        .set{ch_experiment_npooled}

    if (params.run_with_genotype_input) {
        cellsnp_output_dir.combine(ch_experiment_npooled, by: 0)
            .combine(merged_expected_genotypes, by: 0).set{full_vcf}  
    }else{
        log.info "-----running Vireo without genotype input----"
        cellsnp_output_dir.combine(ch_experiment_npooled, by: 0).set{full_vcf}
        full_vcf.map {experiment, cellsnp, npooled -> tuple(experiment, cellsnp, npooled,[],[])}.set{full_vcf}
    }

    full_vcf.filter { experiment, cellsnp, npooled, t,ti -> npooled != '1' }.set{full_vcf2}
    full_vcf.filter { experiment, cellsnp, npooled, t,ti -> npooled == '1' }.set{not_deconvoluted}
            
    VIREO(full_vcf2)  
    // To make things constent we still use donor0 ... donor1 input names, hence the folowing module to change this back to a default approach
    
    REPLACE_GT_DONOR_ID2(VIREO.out.all_required_data)

    vireo_out_sample_donor_vcf = REPLACE_GT_DONOR_ID2.out.infered_vcf
    vireo_out_sample_summary_tsv = REPLACE_GT_DONOR_ID2.out.sample_summary_tsv
    vireo_out_sample__exp_summary_tsv = REPLACE_GT_DONOR_ID2.out.sample__exp_summary_tsv
    vireo_out_sample_donor_ids = REPLACE_GT_DONOR_ID2.out.sample_donor_ids

    // vireo_out_sample_donor_vcf2 = VIREO.out.infered_vcf
    // vireo_out_sample_summary_tsv2 = VIREO.out.sample_summary_tsv
    // vireo_out_sample__exp_summary_tsv2 = VIREO.out.sample__exp_summary_tsv
    // vireo_out_sample_donor_ids2 = VIREO.out.sample_donor_ids

    not_deconvoluted.map{ experiment, donorsvcf, npooled,t,t2 -> tuple(experiment, 'None')}.set{not_deconvoluted2}
    file_cellmetadata = MULTIPLET.out.file__cellmetadata
    scrublet_paths = MULTIPLET.out.scrublet_paths
    split_channel = vireo_out_sample_donor_ids.combine(ch_experiment_filth5, by: 0)
    split_channel2 = not_deconvoluted2.combine(ch_experiment_filth5, by: 0)
    // combining these 2 channels in one
    split_channel3 = split_channel.mix(split_channel2)
    // adding the scrublet paths to the channel.
    split_channel4 = split_channel3.combine(scrublet_paths, by: 0)

    split_channel5 = split_channel4.map{
        val_sample, val_donor_ids_tsv, val_filtered_matrix_h5, path_scrublet ->
        [  val_sample,file(val_donor_ids_tsv),file(val_filtered_matrix_h5),path_scrublet,params.outdir]
    }
    
    SPLIT_DONOR_H5AD(split_channel5)
    if (params.run_with_genotype_input) {
        match_genotypes(vireo_out_sample_donor_vcf,merged_expected_genotypes)
    }

}


workflow REPORT_UPDATE{
    // We use this entry point to update the reports upon running some individual processes that have already completed.
    input_channel = Channel.fromPath(params.input_data_table, followLinks: true, checkIfExists: true)
    CREATE_ARTIFICIAL_BAM_CHANNEL(input_channel)
    bam_split_channel = CREATE_ARTIFICIAL_BAM_CHANNEL.out.ch_experiment_bam_bai_barcodes
    ch_poolid_csv_donor_assignments = CREATE_ARTIFICIAL_BAM_CHANNEL.out.ch_poolid_csv_donor_assignments
    

    // 1) The pihat values were impemented posthoc, hence we are runing this on each of the independent tranches. 
    myFileChannel = Channel.fromPath( "${params.outdir}/deconvolution/vireo/*/GT_donors.vireo.vcf.gz" )
    myFileChannel.map{row -> tuple(row[-2], row)}.set{vireo_out_sample_donor_vcf}
    match_genotypes(vireo_out_sample_donor_vcf)
    match_genotypes.out.out_finish_val.set{o1}
    // updating the Metadata if something new has been fetched,
    // UKBB is sending us samples once a week and sometimes the sample mappings may be present at a later date, hence we update previously run samples accordingly.
    Channel.from([["${params.RUN}","${params.output_dir}"]]).set{update_input_channel}
    // We sometimes aslo chnge apporach in the data fetch and we need to add in some extra metadata
    // metadata_posthoc(update_input_channel)
    replace_donors_posthoc(update_input_channel)
    replace_donors_posthoc.out.dummy_out.set{o2}
    o1.mix(o2).last().set{o3}
    // Once everything is updated we need to make sure that the dataon the website and in the cardinal analysis foder is accurate and up to date, hence we rerun the data_handover scripts.
    // data_handover(params.output_dir,
    //             process_finish_check_channel,
    //             ch_poolid_csv_donor_assignments,
    //             bam_split_channel) 
    
    SUMMARY_STATISTICS_PLOTS(params.output_dir,o3,params.input_data_table)
    TRANSFER(SUMMARY_STATISTICS_PLOTS.out.summary_plots,params.rsync_to_web_file,params.output_dir)
}


/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { SCDECON } from './workflows/yascp'

//
// WORKFLOW: Run main nf-core/yascp analysis pipeline
//
workflow NFCORE_SCDECON {
    SCDECON ()
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_SCDECON ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
