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
        mode            // 'default' or upstream assignment data
        input_channel   // Raw input table
        vcf_input       // Donor VCFs for deconvolution

    main:
        ch_versions = Channel.empty()

        // 1. INITIALIZE ASSIGNMENTS & INPUTS
        if ("${mode}" != 'default') {
            // here we have rerun something upstream - done for freeze1
            assignments_all_pools = mode
        }

        PREPARE_INPUTS(input_channel)

        if (!params.input_data_table.contains('fake_file')) {
            log.info " ---- Genome used: ${params.reference_assembly_fasta_dir} ---"
            input_channel = PREPARE_INPUTS.out.channel_input_data_table
            
            if (params.reference_assembly_fasta_dir == 'https://yascp.cog.sanger.ac.uk/public/10x_reference_assembly') {
                RETRIEVE_RECOURSES()  
                genome1 = RETRIEVE_RECOURSES.out.reference_assembly
            } else {
                genome1 = "${params.reference_assembly_fasta_dir}"
            }

            PREPROCESS_GENOME(genome1)
            genome         = PREPROCESS_GENOME.out.preprocessed_genome
            chanel_cr_outs = PREPARE_INPUTS.out.chanel_cr_outs
            channel_dsb    = PREPARE_INPUTS.out.channel_dsb
            ch_versions    = ch_versions.mix(PREPROCESS_GENOME.out.versions)
        }

        // Initialize placeholders
        vireo_paths                     = Channel.from("$projectDir/assets/fake_file.fq")
        matched_donors                  = Channel.from("$projectDir/assets/fake_file.fq")
        ch_poolid_csv_donor_assignments = Channel.empty()
        bam_split_channel               = Channel.of()
        out_ch = params.outdir
            ? Channel.fromPath(params.outdir, checkIfExists:true)
            : Channel.from("${launchDir}/${params.outdir}")

        // 2. MAIN PROCESSING BRANCH
        if (!params.just_reports) {
            
            if (!params.skip_preprocessing) {
                // CITESEQ and other data modality separation
                SPLIT_CITESEQ_GEX_FILTERED(PREPARE_INPUTS.out.ch_experimentid_paths10x_filtered, 'filtered')
                SPLIT_CITESEQ_GEX(PREPARE_INPUTS.out.ch_experimentid_paths10x_raw, 'raw')

                // Ambient RNA removal (Cellbender vs Cellranger)
                if (params.input == 'cellbender') {
                    if (params.cellbender_with_citeseq) {
                        log.info ' ---- Running cellbender with citeseq ---'
                        ch_experimentid_paths10x_raw = PREPARE_INPUTS.out.ch_experimentid_paths10x_raw
                    } else {
                        ch_experimentid_paths10x_raw = SPLIT_CITESEQ_GEX.out.gex_data
                    }

                    AMBIENT_RNA(ch_experimentid_paths10x_raw, PREPARE_INPUTS.out.ch_experimentid_paths10x_filtered, PREPARE_INPUTS.out.channel__metadata)
                    ch_versions = ch_versions.mix(AMBIENT_RNA.out.cellbender_versions)

                    SPLIT_CITESEQ_GEX_FILTERED_NOCB(AMBIENT_RNA.out.cellbender_path, 'filtered_after_cb')
                    SPLIT_CITESEQ_GEX_NOCB(AMBIENT_RNA.out.cellbender_path_raw, 'raw_after_cb')
                    
                    DECONV_INPUTS(AMBIENT_RNA.out.cellbender_path, PREPARE_INPUTS)

                    channel__file_paths_10x        = DECONV_INPUTS.out.channel__file_paths_10x
                    ch_experiment_bam_bai_barcodes = DECONV_INPUTS.out.ch_experiment_bam_bai_barcodes
                    ch_experiment_filth5           = AMBIENT_RNA.out.cellbender_path
                }
                else if (params.input == 'cellranger') {
                    log.info '--- using cellranger filtered data instead of cellbender (skipping cellbender)---'
                    channel__file_paths_10x        = SPLIT_CITESEQ_GEX_FILTERED.out.gex_data
                    ch_experiment_filth5           = SPLIT_CITESEQ_GEX_FILTERED.out.gex_data
                    ch_experiment_bam_bai_barcodes = PREPARE_INPUTS.out.ch_experiment_bam_bai_barcodes
                }
                else {
                    log.error '--- input mode is not selected - please choose (cellbender | cellranger) ---'
                }

                // HASTAG DECONVOLUTION
                SPLIT_CITESEQ_GEX.out.multiplexing_capture_channel_for_demultiplexing
                    .map { sample_name, path1, path2 ->
                        def directories = path1.findAll { it.isDirectory() }
                        tuple(sample_name, directories, path2)
                    }
                    .set { filtered_multiplexing_capture_channel }

                HASTAG_DEMULTIPLEX(filtered_multiplexing_capture_channel)
                HASTAG_FILE_MERGE(HASTAG_DEMULTIPLEX.out.results.collect(), 'hastag')
                hastag_labels = HASTAG_FILE_MERGE.out.results.ifEmpty { "$projectDir/assets/fake_file1.fq" }

                // MODALITY SELECTION FOR DOWNSTREAM
                if (params.doublets_and_celltypes_on_cellbender_corrected_counts && params.input == 'cellbender') {
                    channel__file_paths_10x_gex = SPLIT_CITESEQ_GEX_FILTERED_NOCB.out.gex_data
                } else {
                    channel__file_paths_10x_gex = SPLIT_CITESEQ_GEX_FILTERED.out.gex_data
                }

                // DOUBLET DETECTION
                if (params.filter_multiplets.run_process) {
                    log.info '---Running doublet assignment ----'
                    MULTIPLET(channel__file_paths_10x_gex, 'yascp_full')
                    doublet_paths = MULTIPLET.out.scrublet_paths
                    sf_mult       = MULTIPLET.out.result_sf.collect()
                    ch_versions   = ch_versions.mix(MULTIPLET.out.doublet_versions)
                } else {
                    log.info '---Skipping doublet assignment ----'
                    sf_mult       = Channel.of()
                    doublet_paths = Channel.from("$projectDir/assets/fake_file.fq")
                }

                DOUBLET_FILE_MERGE(sf_mult, 'doublet')
                doublet_labels = DOUBLET_FILE_MERGE.out.results.ifEmpty { "$projectDir/assets/fake_file2.fq" }

                // CELLTYPE ASSIGNMENT
                if (params.celltype_assignment.run_celltype_assignment) {
                    log.info '---Running celltype assignment ----'
                    CELLTYPE(channel__file_paths_10x_gex, 'yascp_full')
                    celltype_assignments = CELLTYPE.out.celltype_assignments
                    ch_versions          = ch_versions.mix(CELLTYPE.out.celltype_versions)
                } else {
                    log.info '---Skipping celltype assignment ----'
                    celltype_assignments = Channel.from("$projectDir/assets/fake_file.fq")
                }

                // GENETIC DECONVOLUTION
                // When thepreprocessing with cellbender or cellranger is finalised then we can do the deconvolution of samples. This can also be skipped if the samples are not multiplexed.
                // However if the number of individuals is specified as 1 the deconvolution will be skipped anyways, but we will apply scrubblet to remove dublicates.

                if (params.do_deconvolution) {
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

                    ch_versions                     = ch_versions.mix(MAIN_DECONVOLUTION.out.deconvolution_versions)
                    vireo_paths                     = MAIN_DECONVOLUTION.out.vireo_paths2
                    matched_donors                  = MAIN_DECONVOLUTION.out.matched_donors
                    ch_poolid_csv_donor_assignments = MAIN_DECONVOLUTION.out.ch_poolid_csv_donor_assignments
                    bam_split_channel               = MAIN_DECONVOLUTION.out.sample_possorted_bam_vireo_donor_ids
                    assignments_all_pools           = MAIN_DECONVOLUTION.out.assignments_all_pools
                    
                    if (!params.atac) {
                        MERGE_SAMPLES(MAIN_DECONVOLUTION.out.out_h5ad, MAIN_DECONVOLUTION.out.vireo_out_sample__exp_summary_tsv, celltype_assignments, hastag_labels, doublet_labels, 'h5ad')
                    }
                } else {
                    log.info '--- Skipping Deconvolution ---'
                    channel__metadata     = PREPARE_INPUTS.out.channel__metadata
                    MERGE_SAMPLES(channel__file_paths_10x, channel__metadata, celltype_assignments, hastag_labels, doublet_labels, 'barcodes')
                    assignments_all_pools = Channel.from("$projectDir/assets/fake_file1.fq")
                    vireo_paths           = Channel.from("$projectDir/assets/fake_file.fq")
                    matched_donors        = Channel.from("$projectDir/assets/fake_file.fq")
                }
                
                if (!params.atac) {
                    file__anndata_merged = MERGE_SAMPLES.out.file__anndata_merged
                    DUMMY_FILTERED_CHANNEL(file__anndata_merged, params.id_in)
                    file__cells_filtered = DUMMY_FILTERED_CHANNEL.out.anndata_metadata
                }

            } else {
                // SKIP PREPROCESSING: Start from prepared h5ad
                log.info '----Skipping Preprocessing since we already have prepeared h5ad input file----'
                file__anndata_merged  = Channel.from(params.file__anndata_merged)
                assignments_all_pools = Channel.from("$projectDir/assets/fake_file.fq")

                vireo_paths = params.outdir
                    ? Channel.fromPath("${params.outdir}/deconvolution/vireo_raw/*/vireo_*", checkIfExists:true, type: 'dir')
                    : Channel.fromPath("${launchDir}/${params.outdir}/deconvolution/vireo_raw/*/vireo_*", type: 'dir')

                GENOTYPE_MATCHER(vireo_paths.collect())
                matched_donors = GENOTYPE_MATCHER.out.matched_donors

                if ("${mode}" != 'default') {
                    assignments_all_pools = mode
                } else {
                    assignments_all_pools = (params.file__anndata_merged != '') 
                        ? Channel.from(params.gt_match_file) 
                        : Channel.from("$projectDir/assets/fake_file.fq")
                }
                
                if (params.file__cells_filtered == '' && !params.atac) {
                    DUMMY_FILTERED_CHANNEL(file__anndata_merged, params.id_in)
                    file__cells_filtered = DUMMY_FILTERED_CHANNEL.out.anndata_metadata
                } else {
                    file__cells_filtered = Channel.from(params.skip_preprocessing.file__cells_filtered)
                }

                CREATE_ARTIFICIAL_BAM_CHANNEL(input_channel)
                bam_split_channel               = CREATE_ARTIFICIAL_BAM_CHANNEL.out.ch_experiment_bam_bai_barcodes
                ch_poolid_csv_donor_assignments = CREATE_ARTIFICIAL_BAM_CHANNEL.out.ch_poolid_csv_donor_assignments
            }

            // 3. QC & INTEGRATION
            if (!params.skip_qc && !params.atac) {
                gt_outlier_input = (params.gt_match_based_adaptive_qc_exclusion_pattern != '') 
                    ? assignments_all_pools 
                    : Channel.from("$projectDir/assets/fake_file.fq")
                
                QC_AND_INTEGRATION(file__anndata_merged, file__cells_filtered, gt_outlier_input, channel_dsb, vireo_paths, assignments_all_pools, matched_donors, chanel_cr_outs)
                
                ch_versions                  = ch_versions.mix(QC_AND_INTEGRATION.out.qc_versions)
                process_finish_check_channel = QC_AND_INTEGRATION.out.LI
                file__anndata_merged         = QC_AND_INTEGRATION.out.file__anndata_merged
            } else {
                process_finish_check_channel = Channel.of([1, 'dummy'])
            }

        } else {
            // JUST REPORTS branch
            CREATE_ARTIFICIAL_BAM_CHANNEL(input_channel)
            bam_split_channel            = CREATE_ARTIFICIAL_BAM_CHANNEL.out.ch_experiment_bam_bai_barcodes
            process_finish_check_channel = Channel.of([1, 'dummy'])
        }

        // 4. DATA HANDOVER
        if (!params.skip_handover || !params.skip_qc) {
            DATA_HANDOVER(
                out_ch, 
                input_channel,
                process_finish_check_channel,
                ch_poolid_csv_donor_assignments,
                bam_split_channel, 
                genome
            ) 
        }

        // 5. VERSIONS
        softwareVersionsToYAML(ch_versions)
            .collectFile(storeDir: "${params.outdir}", name: 'yascp_software_versions.yml', sort: true, newLine: true)
}