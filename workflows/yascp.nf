/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)



// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist

// Check mandatory parameters
// if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

//
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions' addParams( options: [publish_files : ['tsv':'']] )
include {deconvolution; match_genotypes} from "$projectDir/subworkflows/deconvolution"
include {cellbender} from "$projectDir/subworkflows/cellbender"
include {qc} from "$projectDir/subworkflows/qc"
include {data_handover} from "$projectDir/subworkflows/data_handover"
include { prepare_inputs } from "$projectDir/subworkflows/prepare_inputs"
include {MERGE_SAMPLES} from "$projectDir/modules/nf-core/modules/merge_samples/main"
include {MULTIPLET} from "../modules/nf-core/modules/multiplet/main"
include {dummy_filtered_channel} from "../modules/nf-core/modules/merge_samples/functions"
//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check' addParams( options: [:] )

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow DECONV_INPUTS{
    take:
        cellbender_path
        prepare_inputs
    main:
        cellbender_path
            .splitCsv(header: true, sep: params.input_tables_column_delimiter)
            .map{row->tuple(row.experiment_id, file(cellbender_path+'/../../../../'+row.data_path_10x_format))}
            .set{ch_experiment_filth5} // this channel is used for task 'split_donor_h5ad'
        prepare_inputs.out.ch_experiment_bam_bai_barcodes.map { experiment, bam, bai, barcodes -> tuple(experiment,
                            bam,
                            bai)}.set{pre_ch_experiment_bam_bai_barcodes}
        cellbender_path
            .splitCsv(header: true, sep: params.input_tables_column_delimiter).map{row->tuple(row.experiment_id, file(cellbender_path+'/../../../../'+row.data_path_10x_format+'/barcodes.tsv.gz'))}.set{barcodes}
        channel__file_paths_10x= cellbender_path
            .splitCsv(header: true, sep: params.input_tables_column_delimiter).map{row->tuple(row.experiment_id,
                                                    file(cellbender_path+'/../../../../'+row.data_path_10x_format+'/barcodes.tsv.gz'),
                                                    file(cellbender_path+'/../../../../'+row.data_path_10x_format+'/features.tsv.gz'),
                                                    file(cellbender_path+'/../../../../'+row.data_path_10x_format+'/matrix.mtx.gz'))}
        pre_ch_experiment_bam_bai_barcodes.combine(barcodes, by: 0).set{ch_experiment_bam_bai_barcodes}
    emit:
        channel__file_paths_10x
        ch_experiment_bam_bai_barcodes
        ch_experiment_filth5
}



workflow SCDECON {

    if (!params.skip_preprocessing.value){
        input_channel = Channel.fromPath(params.input_data_table, followLinks: true, checkIfExists: true)
        // prepearing the inputs from a standard 10x dataset folders.
        prepare_inputs(input_channel)
        log.info 'The preprocessing has been already performed, skipping directly to h5ad input'
           // // Removing the background using cellbender which is then used in the deconvolution.
        if (params.input == 'cellbender'){
            log.info ' ---- using cellbender to remove background---'
            cellbender(prepare_inputs.out.ch_experimentid_paths10x_raw,
                prepare_inputs.out.ch_experimentid_paths10x_filtered,prepare_inputs.out.channel__metadata)
            log.info ' ---- Out results - cellbender to remove background---'
            DECONV_INPUTS(cellbender.out.cellbender_path,prepare_inputs)
            channel__file_paths_10x = DECONV_INPUTS.out.channel__file_paths_10x
            ch_experiment_bam_bai_barcodes= DECONV_INPUTS.out.ch_experiment_bam_bai_barcodes
            ch_experiment_filth5= DECONV_INPUTS.out.ch_experiment_filth5
        }else if (params.input == 'existing_cellbender'){
            log.info ' ---- using existing cellbender output for deconvolution---'
            cellbender_location = Channel.from(params.cellbender_location+params.cellbender_filenamePattern+params.cellbender_resolution_to_use+'.tsv')
            cellbender_location.view()
            DECONV_INPUTS(cellbender_location,prepare_inputs)
            channel__file_paths_10x = DECONV_INPUTS.out.channel__file_paths_10x
            ch_experiment_bam_bai_barcodes= DECONV_INPUTS.out.ch_experiment_bam_bai_barcodes
            ch_experiment_filth5= DECONV_INPUTS.out.ch_experiment_filth5
        }
        else if (params.input == 'cellranger'){
            log.info '--- using cellranger filtered data instead of cellbender (skipping cellbender)---'
            ch_experiment_bam_bai_barcodes=prepare_inputs.out.ch_experiment_bam_bai_barcodes
            ch_experiment_filth5=prepare_inputs.out.ch_experiment_filth5
            channel__file_paths_10x=prepare_inputs.out.channel__file_paths_10x
        }
        else{
            log.info '--- input mode is not selected - please choose --- (existing_cellbender| cellbender | cellranger)'
        }

        if (params.do_deconvolution){
            Channel.fromPath(params.reference_genotype_vcf)
            .map { file -> tuple(file, "${file}.tbi")}
            .subscribe { println "TEST_MATCH_GT_VIREO: ${it}" }
            .set { ch_ref_vcf }
            // ch_experiment_filth5.view()
            deconvolution(ch_experiment_bam_bai_barcodes, // activate this to run deconvolution pipeline
                prepare_inputs.out.ch_experiment_npooled,
                ch_experiment_filth5,
                prepare_inputs.out.ch_experiment_donorsvcf_donorslist,channel__file_paths_10x)
                MERGE_SAMPLES(deconvolution.out.out_h5ad,deconvolution.out.vireo_out_sample__exp_summary_tsv,'h5ad')
                match_genotypes(deconvolution.out.vireo_out_sample_donor_vcf, ch_ref_vcf)
        }else{
            channel__metadata = prepare_inputs.out.channel__metadata
            MERGE_SAMPLES(channel__file_paths_10x,channel__metadata,'barcodes')
        }
        // Here add a fundtion to take an extra h5ad and merge it together with the current run.
        file__anndata_merged = MERGE_SAMPLES.out.file__anndata_merged
        file__cells_filtered = MERGE_SAMPLES.out.file__cells_filtered
    }else{
        log.info '''----Skipping Preprocessing since we already have prepeared h5ad input file----'''
        file__anndata_merged = Channel.from(params.skip_preprocessing.file__anndata_merged)
        if (params.skip_preprocessing.file__cells_filtered ==''){
            log.info '''---No cells filtered input'''
            dummy_filtered_channel(file__anndata_merged,params.skip_preprocessing.id_in)
            file__cells_filtered = dummy_filtered_channel.out.anndata_metadata

        }else{
            file__cells_filtered = Channel.from(params.skip_preprocessing.file__cells_filtered)
        }
    }

    qc(file__anndata_merged,file__cells_filtered) //This runs the Clusterring and qc assessments of the datasets.
    // // Performing eQTL mapping.

    data_handover("${workDir}/../${params.output_dir}",qc.out.LI) //This part gathers the plots for the reporting in a Summary folder. If run through gitlab CI it will triger the data transfer to web.


}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
