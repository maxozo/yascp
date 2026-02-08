
include { AZIMUTH; AZIMUTH_ATAC } from "$projectDir/modules/local/azimuth/main"
include { CELLTYPIST } from "$projectDir/modules/local/celltypist/main"
include { SPLIT_BATCH_H5AD } from "$projectDir/modules/local/split_batch_h5ad/main"
include { KERAS_CELLTYPE } from "$projectDir/modules/local/keras_celltype/main"
include { SCPRED } from "$projectDir/modules/local/scpred/main"
include { CONVERT_MTX_TO_H5AD } from "$projectDir/modules/local/convert_h5ad_to_mtx/main"

process CELLTYPE_FILE_MERGE{
    tag "${samplename}"    
    label 'process_high'
    publishDir  path: "${params.outdir}/celltype_assignment/",
            saveAs: {filename ->
                    if (filename.contains("adata.h5ad")) {
                        null
                    } else {
                        filename
                    }
                },
            mode: "${params.copy_mode}",
            overwrite: "true"  

    publishDir  path: "${params.outdir}/handover/merged_h5ad/",
            saveAs: {filename ->
                    if (filename.contains("adata.h5ad")) {
                        filename = "2.celltype_anotated_merged.h5ad"
                    } else {
                        null
                    }
                },
            mode: "${params.copy_mode}",
            overwrite: "true"  


    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
       container "${params.yascp_container_docker}"
    }
    output:
        path("All_Celltype_Assignments.tsv",emit:celltype_assignments)
        path "tranche_celltype_report.tsv"
        path "donor_celltype_report.tsv"

    input:
        path(azimuth_files)
        path(celltypist_paths)
        path(all_other_paths)
    script:
        def merged_files_outpath = workflow.workDir.toString()
        file(merged_files_outpath).mkdirs()
        def azimuth_files_path = "${merged_files_outpath}/azimuth_files.tsv"
        def celltypist_files_path = "${merged_files_outpath}/celltypist_files.tsv"
        def all_other_files_path = "${merged_files_outpath}/other_files.tsv"

        new File(azimuth_files_path).text = azimuth_files.join("\n")
        new File(celltypist_files_path).text = celltypist_paths.join("\n")

        if ("${all_other_paths}" != 'fake_file.fq') {
            new File(all_other_files_path).text = all_other_paths.join("\n")
            other_paths = "--all_other_paths ${all_other_files_path}"
        } else {
            other_paths = ""
        }

        """
            generate_combined_celltype_anotation_file.py --all_azimuth_files ${azimuth_files_path} --all_celltypist_files ${celltypist_files_path} ${other_paths}
        """

}


workflow CELLTYPE{
    
    take:
        file__anndata_merged
        mode
    main:
        // Here we may want to not split it and just pass in an entire h5ad file for annotations.
        // We need a combined h5ad file with all donors to perform further data integrations
        Channel.empty().set { ch_versions }
        if (mode=='yascp_full'){
            file__anndata_merged_post = CONVERT_MTX_TO_H5AD(file__anndata_merged).gex_h5ad
        }else{
            log.info '---Splitting the assignment for each batch---'
            SPLIT_BATCH_H5AD(file__anndata_merged,params.doublet_celltype_split_column)
            SPLIT_BATCH_H5AD.out.sample_file
                .splitCsv(header: true, sep: "\t", by: 1)
                .map{row -> tuple(row.experiment_id, file(row.h5ad_filepath))}.set{file__anndata_merged_post}           
        }
        
        //
        ch_experiment_filth5 = file__anndata_merged_post 
        // az_ch_experiment_filth5 = file__anndata_merged_post


        // Keras celltype assignemt
        if (params.celltype_assignment.run_keras){
            KERAS_CELLTYPE(ch_experiment_filth5,params.celltype_prediction.keras.keras_model,params.celltype_prediction.keras.keras_weights_df) 
            ch_versions = ch_versions.mix(KERAS_CELLTYPE.out.versions)
            all_extra_fields3 = KERAS_CELLTYPE.out.predicted_celltype_labels.collect()
            all_extra_fields = all_extra_fields3.ifEmpty(Channel.from("$projectDir/assets/fake_file.fq"))
        }else{
            all_extra_fields = Channel.from("$projectDir/assets/fake_file.fq")
        }
        
        // AZIMUTH
        if (params.celltype_assignment.run_azimuth){
            if (params.atac){
                // Comented out here because its failing with a missing index file issue caused by software versions - 
                // fixed it by installing a previous version of Seurat and Signac https://github.com/satijalab/azimuth/issues/211. But the container has to be updated to work.
                // AZIMUTH_ATAC(file__anndata_merged,params.mapping_file,Channel.fromList( params.azimuth.celltype_atac_refsets))
                az_out = Channel.from("$projectDir/assets/fake_file1.fq")
            }else{
                AZIMUTH(file__anndata_merged,params.mapping_file,Channel.fromList( params.azimuth.celltype_refsets))
                ch_versions = ch_versions.mix(AZIMUTH.out.versions)
                az_out = AZIMUTH.out.predicted_celltype_labels
                    .ifEmpty { "$projectDir/assets/fake_file1.fq" }
            }

        }else{
            az_out = Channel.from("$projectDir/assets/fake_file1.fq")
            az_out = az_out.ifEmpty(Channel.from("$projectDir/assets/fake_file1.fq"))
        }
        
        // CELLTYPIST
        if (params.celltype_assignment.run_celltypist){
            Channel.fromList(params.celltypist.models)
                .set{ch_celltypist_models}
            CELLTYPIST(ch_experiment_filth5.combine(ch_celltypist_models))
            ch_versions = ch_versions.mix(CELLTYPIST.out.versions)
            ct_out = CELLTYPIST.out.predicted_labels
                .ifEmpty { "$projectDir/assets/fake_file2.fq" }
        }else{
            ct_out = Channel.from("$projectDir/assets/fake_file2.fq")
        }

        // // SCPRED
        if (params.celltype_assignment.run_scpred){
            SCPRED(ch_experiment_filth5,params.scpred.reference)
            ch_versions = ch_versions.mix(SCPRED.out.versions)
            sc_out2 = SCPRED.out.predicted_celltype_labels.collect()
            sc_out = sc_out2.ifEmpty(Channel.of())
        }else{
            sc_out = Channel.of()
        }        
        all_extra_fields2 = all_extra_fields.mix(sc_out)
        CELLTYPE_FILE_MERGE(az_out.collect().unique(),ct_out.collect().unique(),all_extra_fields2.collect().unique()) 
        celltype_assignments=CELLTYPE_FILE_MERGE.out.celltype_assignments
    emit:
        celltype_assignments
        celltype_versions = ch_versions

}