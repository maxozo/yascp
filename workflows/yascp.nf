/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/
include { GET_SOFTWARE_VERSIONS } from "$projectDir/modules/local/get_software_versions" addParams( options: [publish_files : ['tsv':'']] )
include { main_deconvolution } from "$projectDir/subworkflows/main_deconvolution"
include {ambient_RNA} from "$projectDir/subworkflows/ambient_RNA"
include {qc} from "$projectDir/subworkflows/qc"
include {celltype} from "$projectDir/subworkflows/celltype"
include {data_handover} from "$projectDir/subworkflows/data_handover"
include { prepare_inputs } from "$projectDir/subworkflows/prepare_inputs"
include { DECONV_INPUTS } from "$projectDir/subworkflows/prepare_inputs"
include { CREATE_ARTIFICIAL_BAM_CHANNEL } from "$projectDir/modules/local/create_artificial_bam_channel/main"
include {MERGE_SAMPLES} from "$projectDir/modules/nf-core/modules/merge_samples/main"
include {dummy_filtered_channel} from "$projectDir/modules/nf-core/modules/merge_samples/functions"
include {MULTIPLET} from "$projectDir/subworkflows/doublet_detection"
include { SPLIT_CITESEQ_GEX; SPLIT_CITESEQ_GEX as SPLIT_CITESEQ_GEX_FILTERED;SPLIT_CITESEQ_GEX as SPLIT_CITESEQ_GEX_FILTERED_NOCB;SPLIT_CITESEQ_GEX as SPLIT_CITESEQ_GEX_NOCB; SPLIT_CITESEQ_GEX as PREPOCESS_FILES; HASTAG_DEMULTIPLEX } from '../modules/nf-core/modules/citeseq/main'
include { GENOTYPE_MATCHER } from "$projectDir/modules/nf-core/modules/vireo/main"
include { RETRIEVE_RECOURSES } from "$projectDir/subworkflows/local/retrieve_recourses"
include { PREPROCESS_GENOME } from "$projectDir/modules/nf-core/modules/subset_bam_per_barcodes_and_variants/main"
/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// This is the main workflow which consists of:
//  1) Ambient RNA removal using cellbender - this is a lenghty process as GPU is required. ../subworkflows/ambient_RNA.nf
//  2) Deconvolution and GT match (if genotypes provided) ../subworkflows/main_deconvolution.nf
//  3) Celltype assignment  ../subworkflows/celltype.nf
//  5) Data handover preparation  ../subworkflows/data_handover.nf


workflow YASCP {
    take:
        mode
        input_channel
        input_channel
        vcf_input
    main:
        if("${mode}"!='default'){
            // here we have rerun something upstream - done for freeze1
            assignments_all_pools = mode
        }
        prepare_inputs(input_channel)
        if (!params.input_data_table.contains('fake_file')){
            input_channel = prepare_inputs.out.channel_input_data_table
            if (params.reference_assembly_fasta_dir=='https://yascp.cog.sanger.ac.uk/public/10x_reference_assembly'){
                RETRIEVE_RECOURSES()  
                genome1 = RETRIEVE_RECOURSES.out.reference_assembly
            }else{
                genome1 = "${params.reference_assembly_fasta_dir}"
            }
            genome = PREPROCESS_GENOME(genome1)
                
            chanel_cr_outs = prepare_inputs.out.chanel_cr_outs
            channel_dsb = prepare_inputs.out.channel_dsb
        }
        vireo_paths = Channel.from("$projectDir/assets/fake_file.fq")
        matched_donors = Channel.from("$projectDir/assets/fake_file.fq")

        ch_poolid_csv_donor_assignments = Channel.empty()
        bam_split_channel = Channel.of()
        out_ch = params.outdir
            ? Channel.fromPath(params.outdir, checkIfExists:true)
            : Channel.from("${launchDir}/${params.outdir}")
                 
        if(!params.just_reports){
            // sometimes we just want to rerun report generation as a result of alterations, hence if we set params.just_reports =True pipeline will use the results directory and generate a new reports.
            if (!params.skip_preprocessing){
                // The input table should contain the folowing columns - experiment_id	n_pooled	donor_vcf_ids	data_path_10x_format
                // prepearing the inputs from a standard 10x dataset folders.

                log.info 'The preprocessing has been already performed, skipping directly to h5ad input'
                // // Removing the background using cellbender which is then used in the deconvolution.

                // CITESEQ and other data modality seperation
                // Split citeseq if available


                // If citeseq data is present in the 10x mtx then we strip it before the ambient rna correction.
                SPLIT_CITESEQ_GEX_FILTERED(prepare_inputs.out.ch_experimentid_paths10x_filtered,'filterd')
                SPLIT_CITESEQ_GEX( prepare_inputs.out.ch_experimentid_paths10x_raw,'raw')

                // Either run ambient RNA removal with cellbender or use cellranger filtered reads (cellbender|cellranger)
                if (params.input == 'cellbender'){
                    if (params.cellbender_with_citeseq){                          // here we run CB with citeseq
                        log.info ' ---- Running cellbender with citeseq ---'
                        ch_experimentid_paths10x_raw = prepare_inputs.out.ch_experimentid_paths10x_raw
                    }else{
                        ch_experimentid_paths10x_raw = SPLIT_CITESEQ_GEX.out.gex_data
                    }
                    // Here we either run ambient RNA removal with citeseq counts or without.
                    ambient_RNA( ch_experimentid_paths10x_raw,
                        prepare_inputs.out.ch_experimentid_paths10x_filtered,prepare_inputs.out.channel__metadata)

                    // Now we convert the CB processed files to h5ad files and split the modalities if they were left in
                    SPLIT_CITESEQ_GEX_FILTERED_NOCB(ambient_RNA.out.cellbender_path,'filterd_after_cb')
                    SPLIT_CITESEQ_GEX_NOCB( ambient_RNA.out.cellbender_path_raw,'raw_after_cb')
                    
                    DECONV_INPUTS(ambient_RNA.out.cellbender_path,prepare_inputs)

                    channel__file_paths_10x = DECONV_INPUTS.out.channel__file_paths_10x
                    ch_experiment_bam_bai_barcodes= DECONV_INPUTS.out.ch_experiment_bam_bai_barcodes
                    ch_experiment_filth5= ambient_RNA.out.cellbender_path

                }
                else if (params.input == 'cellranger'){
                    // This is where we skip the cellbender and use the cellranger filtered datasets.
                    log.info '--- using cellranger filtered data instead of cellbender (skipping cellbender)---'
                    channel__file_paths_10x = SPLIT_CITESEQ_GEX_FILTERED.out.gex_data
                    ch_experiment_filth5 = SPLIT_CITESEQ_GEX_FILTERED.out.gex_data
                    ch_experiment_bam_bai_barcodes=prepare_inputs.out.ch_experiment_bam_bai_barcodes
                    
                }
                else{
                    log.info '--- input mode is not selected - please choose --- (existing_cellbender cellranger)'
                }


                // If we have multiplexing capture file then we proceed with hastag deconvolution
                SPLIT_CITESEQ_GEX.out.multiplexing_capture_channel_for_demultiplexing
                    .map { sample_name, path1, path2 ->
                        def directories = path1.findAll { it.isDirectory() }
                        tuple(sample_name, directories, path2)
                    }
                    .set { filtered_multiplexing_capture_channel }
                HASTAG_DEMULTIPLEX(filtered_multiplexing_capture_channel)
                hastag_labels = HASTAG_DEMULTIPLEX.out.results
                
                if (params.doublets_and_celltypes_on_cellbender_corrected_counts && params.input == 'cellbender'){
                    channel__file_paths_10x_gex = SPLIT_CITESEQ_GEX_FILTERED_NOCB.out.channel__file_paths_10x
                }
                else{
                    channel__file_paths_10x_gex = SPLIT_CITESEQ_GEX_FILTERED.out.channel__file_paths_10x
                }

                // ###################################
                // ###################################
                // Step: DOUBLET DETECTION
                // Curently contains only Scrublet, but we are also adding DoubletDetect
                // ###################################
                // ###################################
                if (params.filter_multiplets.run_process){
                    MULTIPLET(
                        channel__file_paths_10x_gex,'yascp_full'
                    )
                    doublet_paths = MULTIPLET.out.scrublet_paths
                }else{
                    doublet_paths = Channel.from("$projectDir/assets/fake_file.fq")
                }


                if (params.celltype_assignment.run_celltype_assignment){
                    celltype(channel__file_paths_10x_gex,'yascp_full')
                    celltype_assignments=celltype.out.celltype_assignments
                }else{
                    celltype_assignments = Channel.from("$projectDir/assets/fake_file.fq")
                }

                // ###################################
                // ################################### Readme
                // Step2. DECONVOLUTION
                // When thepreprocessing with cellbender or cellranger is finalised then we can do the deconvolution of samples. This can also be skipped if the samples are not multiplexed.
                // However if the number of individuals is specified as 1 the deconvolution withh be skipped anyways, but we will apply scrubblet to remove dublicates.
                // Suggestion is to still run deconvolution so that dublicates are removed.
                // ###################################
                // ###################################
                
                if (params.do_deconvolution){
                    main_deconvolution(ch_experiment_bam_bai_barcodes, // activate this to run deconvolution pipeline
                        prepare_inputs.out.ch_experiment_npooled,
                        ch_experiment_filth5,
                        prepare_inputs.out.ch_experiment_donorsvcf_donorslist,
                        doublet_paths,
                        vcf_input,
                        genome)
                    vireo_paths = main_deconvolution.out.vireo_paths2
                    matched_donors = main_deconvolution.out.matched_donors
                    ch_poolid_csv_donor_assignments = main_deconvolution.out.ch_poolid_csv_donor_assignments
                    bam_split_channel = main_deconvolution.out.sample_possorted_bam_vireo_donor_ids
                    assignments_all_pools = main_deconvolution.out.assignments_all_pools

                    if (!params.skip_merge){
                        log.info '--- merging samples'
                        MERGE_SAMPLES(main_deconvolution.out.out_h5ad,main_deconvolution.out.vireo_out_sample__exp_summary_tsv,celltype_assignments,'h5ad')
                    }else{
                        file__anndata_merged = main_deconvolution.out.out_h5ad
                        dummy_filtered_channel(file__anndata_merged,params.id_in)
                        file__cells_filtered = dummy_filtered_channel.out.anndata_metadata
                    }
                }else{
                    channel__metadata = prepare_inputs.out.channel__metadata
                    if (!params.skip_merge){
                        MERGE_SAMPLES(channel__file_paths_10x,channel__metadata,celltype_assignments,'barcodes')
                    }
                    assignments_all_pools = Channel.from("$projectDir/assets/fake_file.fq")
                    vireo_paths = Channel.from("$projectDir/assets/fake_file.fq")
                    matched_donors = Channel.from("$projectDir/assets/fake_file.fq")
                }
                
                if (!params.skip_merge){
                    file__anndata_merged = MERGE_SAMPLES.out.file__anndata_merged
                    dummy_filtered_channel(file__anndata_merged,params.id_in)
                    file__cells_filtered = dummy_filtered_channel.out.anndata_metadata
                }
            }else{
                // This option skips all the deconvolution and and takes a preprocessed yascp h5ad file to run the downstream clustering and celltype annotation.
                log.info '''----Skipping Preprocessing since we already have prepeared h5ad input file----'''
                file__anndata_merged = Channel.from(params.file__anndata_merged)
                assignments_all_pools = Channel.from("$projectDir/assets/fake_file.fq")

                if (params.citeseq){

                    vireo_paths = params.outdir
                        ? Channel.fromPath("${params.outdir}/deconvolution/vireo_raw/*/vireo_*", checkIfExists:true, type: 'dir')
                        : Channel.fromPath("${launchDir}/${params.outdir}/deconvolution/vireo_raw/*/vireo_*", type: 'dir')

                    GENOTYPE_MATCHER(vireo_paths.collect())
                    matched_donors = GENOTYPE_MATCHER.out.matched_donors
                }else{
                    vireo_paths = Channel.from("$projectDir/assets/fake_file.fq")
                    matched_donors = Channel.from("$projectDir/assets/fake_file.fq")
                }
                
                if("${mode}"!='default'){
                    // Here we have rerun GT matching upstream - done for freeze1
                    assignments_all_pools = mode
                }else{
                    if (params.file__anndata_merged !=''){
                        assignments_all_pools = Channel.from(params.gt_match_file)
                    }else{
                        assignments_all_pools = Channel.from("$projectDir/assets/fake_file.fq")
                    }
                }
                
                if (params.file__cells_filtered ==''){
                    log.info '''--- No cells filtered input ----'''
                    dummy_filtered_channel(file__anndata_merged,params.id_in)
                    file__cells_filtered = dummy_filtered_channel.out.anndata_metadata
                }else{
                    file__cells_filtered = Channel.from(params.skip_preprocessing.file__cells_filtered)
                }
                CREATE_ARTIFICIAL_BAM_CHANNEL(input_channel)
                bam_split_channel = CREATE_ARTIFICIAL_BAM_CHANNEL.out.ch_experiment_bam_bai_barcodes
                ch_poolid_csv_donor_assignments = CREATE_ARTIFICIAL_BAM_CHANNEL.out.ch_poolid_csv_donor_assignments

            }


            // ###################################
            // ################################### Readme
            // QC METRICS, CLUSTERING
            // After background removal and demultiplexing we perform qc metrics and clustering of the processed cells.
            // This step of the pipeline also performs celltype assignments and removes cells that fail adaptive filtering.
            // ###################################
            // ###################################

            if (!params.skip_qc){

                if(params.gt_match_based_adaptive_qc_exclusion_pattern !=''){
                    gt_outlier_input = assignments_all_pools
                }else{
                    gt_outlier_input = Channel.from("$projectDir/assets/fake_file.fq")
                }

                qc(file__anndata_merged,file__cells_filtered,gt_outlier_input,channel_dsb,vireo_paths,assignments_all_pools,matched_donors,chanel_cr_outs) //This runs the Clusterring and qc assessments of the datasets.
                process_finish_check_channel = qc.out.LI
                file__anndata_merged = qc.out.file__anndata_merged
            }else{
                // if we are not running qc step we need to account for an dummy channel. 
                process_finish_check_channel = Channel.of([1, 'dummy'])
            }

        }else{
            // since for the downstreem preocess we do a bam split, and this is generated as part of a main_deconvolution step, we have to generate this input artificially here based on the results directory and fech location.
            CREATE_ARTIFICIAL_BAM_CHANNEL(input_channel)
            bam_split_channel = CREATE_ARTIFICIAL_BAM_CHANNEL.out.ch_experiment_bam_bai_barcodes
            process_finish_check_channel = Channel.of([1, 'dummy'])

        }

        // ###################################
        // ################################### Readme
        // DATA HANDOVER, REPORTS, DATA ENCRYPTION, DONOR H5AD, BAM SPLIT
        // 
        // ###################################
        // ###################################

        if (!params.skip_handover){

            data_handover(out_ch,input_channel,
                            process_finish_check_channel,
                            ch_poolid_csv_donor_assignments,
                            bam_split_channel,genome) 
        }
                        
                        
}

/*
========================================================================================
    THE END
========================================================================================
*/
