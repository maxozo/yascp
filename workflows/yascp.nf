/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

include { MAIN_DECONVOLUTION } from "$projectDir/subworkflows/main_deconvolution"
include {AMBIENT_RNA} from "$projectDir/subworkflows/ambient_RNA"
include {QC_AND_INTEGRATION} from "$projectDir/subworkflows/qc_and_integration"
include {CELLTYPE} from "$projectDir/subworkflows/celltype"
include {DATA_HANDOVER} from "$projectDir/subworkflows/data_handover"
include { PREPARE_INPUTS } from "$projectDir/subworkflows/prepare_inputs"
include { DECONV_INPUTS } from "$projectDir/subworkflows/prepare_inputs"
include { CREATE_ARTIFICIAL_BAM_CHANNEL } from "$projectDir/modules/local/create_artificial_bam_channel/main"
include {MERGE_SAMPLES; HASTAG_FILE_MERGE; HASTAG_FILE_MERGE as DOUBLET_FILE_MERGE} from "$projectDir/modules/local/merge_samples/main"
include {DUMMY_FILTERED_CHANNEL} from "$projectDir/modules/local/merge_samples/functions"
include {MULTIPLET} from "$projectDir/subworkflows/doublet_detection"
include { SPLIT_CITESEQ_GEX; SPLIT_CITESEQ_GEX as SPLIT_CITESEQ_GEX_FILTERED;SPLIT_CITESEQ_GEX as SPLIT_CITESEQ_GEX_FILTERED_NOCB;SPLIT_CITESEQ_GEX as SPLIT_CITESEQ_GEX_NOCB; SPLIT_CITESEQ_GEX as PREPOCESS_FILES; HASTAG_DEMULTIPLEX } from '../modules/local/citeseq/main'
include { GENOTYPE_MATCHER } from "$projectDir/modules/local/vireo/main"
include { RETRIEVE_RECOURSES } from "$projectDir/modules/local/retrieve_resources/retrieve_resources"
include { PREPROCESS_GENOME } from "$projectDir/modules/local/subset_bam_per_barcodes_and_variants/main"
include { softwareVersionsToYAML} from "$projectDir/subworkflows/utils"
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
        vcf_input
    main:
        Channel.empty().set { ch_versions }
        if("${mode}"!='default'){
            // here we have rerun something upstream - done for freeze1
            assignments_all_pools = mode
        }
        PREPARE_INPUTS(input_channel)
        if (!params.input_data_table.contains('fake_file')){
            log.info " ---- Genome used: ${params.reference_assembly_fasta_dir} ---"
            input_channel = PREPARE_INPUTS.out.channel_input_data_table
            if (params.reference_assembly_fasta_dir=='https://yascp.cog.sanger.ac.uk/public/10x_reference_assembly'){
                RETRIEVE_RECOURSES()  
                genome1 = RETRIEVE_RECOURSES.out.reference_assembly
            }else{
                genome1 = "${params.reference_assembly_fasta_dir}"
            }

            PREPROCESS_GENOME(genome1)
            genome = PREPROCESS_GENOME.out.preprocessed_genome
            chanel_cr_outs = PREPARE_INPUTS.out.chanel_cr_outs
            channel_dsb = PREPARE_INPUTS.out.channel_dsb
            ch_versions = ch_versions.mix(PREPROCESS_GENOME.out.versions)
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
                // CITESEQ and other data modality seperation
                // Split citeseq if available
                // If citeseq data is present in the 10x mtx then we strip it before the ambient rna correction.
                SPLIT_CITESEQ_GEX_FILTERED(PREPARE_INPUTS.out.ch_experimentid_paths10x_filtered,'filtered')
                SPLIT_CITESEQ_GEX( PREPARE_INPUTS.out.ch_experimentid_paths10x_raw,'raw')

                // Either run ambient RNA removal with cellbender or use cellranger filtered reads (cellbender|cellranger)
                if (params.input == 'cellbender'){
                    if (params.cellbender_with_citeseq){                          // here we run CB with citeseq
                        log.info ' ---- Running cellbender with citeseq ---'
                        ch_experimentid_paths10x_raw = PREPARE_INPUTS.out.ch_experimentid_paths10x_raw
                    }else{
                        ch_experimentid_paths10x_raw = SPLIT_CITESEQ_GEX.out.gex_data
                    }
                    // Here we either run ambient RNA removal with citeseq counts or without.
                    AMBIENT_RNA( ch_experimentid_paths10x_raw,
                        PREPARE_INPUTS.out.ch_experimentid_paths10x_filtered,PREPARE_INPUTS.out.channel__metadata)
                    ch_versions = ch_versions.mix(AMBIENT_RNA.out.cellbender_versions)

                    // Now we convert the CB processed files to h5ad files and split the modalities if they were left in
                    SPLIT_CITESEQ_GEX_FILTERED_NOCB(AMBIENT_RNA.out.cellbender_path,'filtered_after_cb')
                    SPLIT_CITESEQ_GEX_NOCB( AMBIENT_RNA.out.cellbender_path_raw,'raw_after_cb')
                    
                    DECONV_INPUTS(AMBIENT_RNA.out.cellbender_path,PREPARE_INPUTS)

                    channel__file_paths_10x = DECONV_INPUTS.out.channel__file_paths_10x
                    ch_experiment_bam_bai_barcodes= DECONV_INPUTS.out.ch_experiment_bam_bai_barcodes
                    ch_experiment_filth5= AMBIENT_RNA.out.cellbender_path

                }
                else if (params.input == 'cellranger'){
                    // This is where we skip the cellbender and use the cellranger filtered datasets.
                    log.info '--- using cellranger filtered data instead of cellbender (skipping cellbender)---'
                    channel__file_paths_10x = SPLIT_CITESEQ_GEX_FILTERED.out.gex_data
                    ch_experiment_filth5 = SPLIT_CITESEQ_GEX_FILTERED.out.gex_data
                    ch_experiment_bam_bai_barcodes=PREPARE_INPUTS.out.ch_experiment_bam_bai_barcodes
                }
                else{
                    log.info '--- input mode is not selected - please choose --- (existing_cellbender cellranger)'X
                }

                // If we have multiplexing capture file then we proceed with hastag deconvolution
                SPLIT_CITESEQ_GEX.out.multiplexing_capture_channel_for_demultiplexing
                    .map { sample_name, path1, path2 ->
                        def directories = path1.findAll { it.isDirectory() }
                        tuple(sample_name, directories, path2)
                    }
                    .set { filtered_multiplexing_capture_channel }

                HASTAG_DEMULTIPLEX(filtered_multiplexing_capture_channel)
                HASTAG_FILE_MERGE(HASTAG_DEMULTIPLEX.out.results.collect(),'hastag')
                hastag_labels = HASTAG_FILE_MERGE.out.results
                    .ifEmpty { "$projectDir/assets/fake_file1.fq" }
                

                if (params.doublets_and_celltypes_on_cellbender_corrected_counts && params.input == 'cellbender'){
                    channel__file_paths_10x_gex = SPLIT_CITESEQ_GEX_FILTERED_NOCB.out.gex_data
                }
                else{
                    channel__file_paths_10x_gex = SPLIT_CITESEQ_GEX_FILTERED.out.gex_data
                }
                // channel__file_paths_10x_gex.subscribe { println "channel__file_paths_10x_gex: $it" }
                // // ###################################
                // ###################################
                // Step: DOUBLET DETECTION
                // Curently contains only Scrublet, but we are also adding DoubletDetect
                // ###################################
                // ###################################
                if (params.filter_multiplets.run_process){
                    log.info '---Running doublet assignment ----'
                    MULTIPLET(channel__file_paths_10x_gex,'yascp_full')
                    doublet_paths = MULTIPLET.out.scrublet_paths
                    sf_mult_p = MULTIPLET.out.result_sf
                    sf_mult=sf_mult_p.collect()
                    ch_versions = ch_versions.mix(MULTIPLET.out.doublet_versions)
                }else{
                    log.info '---Skipping celltype assignment ----'
                    sf_mult = Channel.of()
                    doublet_paths = Channel.from("$projectDir/assets/fake_file.fq")
                }

                DOUBLET_FILE_MERGE(sf_mult,'doublet')
                doublet_labels = DOUBLET_FILE_MERGE.out.results
                    .ifEmpty { "$projectDir/assets/fake_file2.fq" }

                if (params.celltype_assignment.run_celltype_assignment){
                    log.info '---Running celltype assignment ----'
                    CELLTYPE(channel__file_paths_10x_gex,'yascp_full')
                    celltype_assignments=CELLTYPE.out.celltype_assignments
                    ch_versions = ch_versions.mix(CELLTYPE.out.celltype_versions)
                }else{
                    log.info '---Skipping celltype assignment ----'
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
                    log.info '--- Performing Deconvolution ---'
                    MAIN_DECONVOLUTION(
                        ch_experiment_bam_bai_barcodes,
                        PREPARE_INPUTS.out.ch_experiment_npooled,
                        ch_experiment_filth5,
                        PREPARE_INPUTS.out.ch_experiment_donorsvcf_donorslist,
                        doublet_paths,
                        vcf_input,
                        genome
                    )

                    ch_versions = ch_versions.mix(MAIN_DECONVOLUTION.out.deconvolution_versions)
                    vireo_paths = MAIN_DECONVOLUTION.out.vireo_paths2
                    matched_donors = MAIN_DECONVOLUTION.out.matched_donors
                    ch_poolid_csv_donor_assignments = MAIN_DECONVOLUTION.out.ch_poolid_csv_donor_assignments
                    bam_split_channel = MAIN_DECONVOLUTION.out.sample_possorted_bam_vireo_donor_ids
                    assignments_all_pools = MAIN_DECONVOLUTION.out.assignments_all_pools
                    if (!params.atac){
                        MERGE_SAMPLES(MAIN_DECONVOLUTION.out.out_h5ad,MAIN_DECONVOLUTION.out.vireo_out_sample__exp_summary_tsv,celltype_assignments,hastag_labels,doublet_labels,'h5ad')
                    }
                    
                }else{
                    log.info '--- Skipping Deconvolution ---'
                    channel__metadata = PREPARE_INPUTS.out.channel__metadata
                    MERGE_SAMPLES(channel__file_paths_10x,channel__metadata,celltype_assignments,hastag_labels,doublet_labels,'barcodes')
                    assignments_all_pools = Channel.from("$projectDir/assets/fake_file1.fq")
                    vireo_paths = Channel.from("$projectDir/assets/fake_file.fq")
                    matched_donors = Channel.from("$projectDir/assets/fake_file.fq")
                }
                
                if (!params.atac){
                    file__anndata_merged = MERGE_SAMPLES.out.file__anndata_merged
                    DUMMY_FILTERED_CHANNEL(file__anndata_merged,params.id_in)
                    file__cells_filtered = DUMMY_FILTERED_CHANNEL.out.anndata_metadata
                }

            }else{
                // This option skips all the deconvolution and and takes a preprocessed yascp h5ad file to run the downstream clustering and celltype annotation.
                log.info '''----Skipping Preprocessing since we already have prepeared h5ad input file----'''
                file__anndata_merged = Channel.from(params.file__anndata_merged)
                assignments_all_pools = Channel.from("$projectDir/assets/fake_file.fq")

                vireo_paths = params.outdir
                    ? Channel.fromPath("${params.outdir}/deconvolution/vireo_raw/*/vireo_*", checkIfExists:true, type: 'dir')
                    : Channel.fromPath("${launchDir}/${params.outdir}/deconvolution/vireo_raw/*/vireo_*", type: 'dir')

                GENOTYPE_MATCHER(vireo_paths.collect())
                matched_donors = GENOTYPE_MATCHER.out.matched_donors

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
                
                if (params.file__cells_filtered =='' && !params.atac){
                    log.info '''--- No cells filtered input ----'''
                    DUMMY_FILTERED_CHANNEL(file__anndata_merged,params.id_in)
                    file__cells_filtered = DUMMY_FILTERED_CHANNEL.out.anndata_metadata
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

            if (!params.skip_qc && !params.atac){
                if(params.gt_match_based_adaptive_qc_exclusion_pattern !=''){
                    gt_outlier_input = assignments_all_pools
                }else{
                    gt_outlier_input = Channel.from("$projectDir/assets/fake_file.fq")
                }
                QC_AND_INTEGRATION(file__anndata_merged,file__cells_filtered,gt_outlier_input,channel_dsb,vireo_paths,assignments_all_pools,matched_donors,chanel_cr_outs) //This runs the Clusterring and qc assessments of the datasets.
                ch_versions = ch_versions.mix(QC_AND_INTEGRATION.out.qc_versions)
                process_finish_check_channel = QC_AND_INTEGRATION.out.LI
                file__anndata_merged = QC_AND_INTEGRATION.out.file__anndata_merged
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
        // ###################################
        // ###################################

        if (!params.skip_handover || !params.skip_qc ){
            DATA_HANDOVER(out_ch,input_channel,
                            process_finish_check_channel,
                            ch_poolid_csv_donor_assignments,
                            bam_split_channel,genome) 
        }

/*
========================================================================================
   versions
========================================================================================
*/
        version_yaml = Channel.empty()
        version_yaml = softwareVersionsToYAML(ch_versions)
            .collectFile(storeDir: "${params.outdir}", name: 'yascp_software_versions.yml', sort: true, newLine: true)
}

/*
========================================================================================
    THE END
========================================================================================
*/
