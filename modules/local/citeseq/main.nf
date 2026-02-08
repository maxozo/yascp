
process SPLIT_CITESEQ_GEX {
    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }

    publishDir  path: "${params.outdir}/preprocessing/data_modalities_split/${mode}/${sample_name}",
    mode: "${params.copy_mode}",
    overwrite: "true"

    input:
        tuple val(sample_name),path(cellranger_raw) 
        val(mode)

    output:
        tuple val(sample_name), path("${sample_name}__*"), path("*__Multiplexing_Capture.tsv"), emit: multiplexing_capture_channel_for_demultiplexing  optional true
        tuple val(sample_name), path("${sample_name}__Gene_Expression"), emit:gex_data
        tuple val(sample_name), path("antibody-${sample_name}.h5ad"), emit: ab_data2 optional true
        tuple val(sample_name), path("Gene_Expression-${sample_name}.h5ad"), emit: gex_h5ad optional true
        path("Gene_Expression-${sample_name}.h5ad"), emit: gex_h5ad_2 optional true
        path("*.tsv"), emit: quants_data optional true
        tuple val(sample_name), path("${sample_name}__*"), emit: ab_data optional true
        tuple val(sample_name), path("${sample_name}__Gene_Expression/barcodes.tsv.gz"), path("${sample_name}__Gene_Expression/features.tsv.gz"), path("${sample_name}__Gene_Expression/matrix.mtx.gz"), emit: channel__file_paths_10x
 
    script:

        """
            matrix_file="${cellranger_raw}/matrix.mtx"
            compressed_file="${cellranger_raw}/matrix.mtx.gz"

            # Check if the compressed file exists
            if [ ! -f "\$compressed_file" ]; then
                echo "\$compressed_file does not exist. Compressing \$matrix_file..."
            
                # Compress the file without deleting the original
                gzip -c "\$matrix_file" > "\$compressed_file"
            
                echo "Compression complete. \$matrix_file has been compressed to \$compressed_file"
            else
                echo "\$compressed_file already exists. No action needed."
            fi

            matrix_file="${cellranger_raw}/barcodes.tsv"
            compressed_file="${cellranger_raw}/barcodes.tsv.gz"
            # Check if the compressed file exists
            if [ ! -f "\$compressed_file" ]; then
                echo "\$compressed_file does not exist. Compressing \$matrix_file..."
            
                # Compress the file without deleting the original
                gzip -c "\$matrix_file" > "\$compressed_file"
            
                echo "Compression complete. \$matrix_file has been compressed to \$compressed_file"
            else
                echo "\$compressed_file already exists. No action needed."
            fi

            features_file="${cellranger_raw}/features.tsv.gz"
            peaks_file="${cellranger_raw}/peaks.bed"

            # Check if the features.tsv.gz file exists
            if [ ! -f "\$features_file" ]; then
            echo "\$features_file does not exist. Creating it from \$peaks_file..."

            # Create the features.tsv file from the peaks.bed file
            awk 'BEGIN{OFS="\t"} {print \$1 ":" \$2 "-" \$3, \$1 ":" \$2 "-" \$3, "Gene Expression"}' "\$peaks_file" | gzip > "\$features_file"

            echo "Creation of \$features_file is complete."
            else
            echo "\$features_file already exists. No action needed."
            fi

            strip_citeseq.py --raw_data ${cellranger_raw} -o ${sample_name} -ha ${params.citeseq_config.hastag_multiplexing_capture_labels}
        """
}


process HASTAG_DEMULTIPLEX {
    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }

    publishDir  path: "${params.outdir}/deconvolution/hastag_demultiplex/${sample_name}", mode: "${params.copy_mode}",
      overwrite: "true"

    input:
        tuple val(sample_name), path(paths), path(multiplexing_capture)
    output:
        path("${sample_name}__hastag_demux_results.tsv"), emit: results
    script:
        """
            hastag_demultiplex.R
            ln -s hastag_demux_results.tsv ${sample_name}__hastag_demux_results.tsv
        """
}

process DSB_INTEGRATE{

    label 'process_medium'
    tag "${sample_name}"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }

    publishDir  path: "${params.outdir}/citeseq/all_data_integrated", mode: "${params.copy_mode}",
      overwrite: "true"
    
    output:
        // path("out"), emit: outdir
        path("*all_samples_integrated.RDS"), emit: tmp_rds_file
        path("figures__*/*"), emit: figs
        


    input:
        path(tmp_rsd)
        each vars_to_regress
        val(k_anchor)
        val(dims)
        val(ndim_sct)
        val(ndim_citeBgRemoved)
        val(ndim_cite_integrated)
    script:

        if (vars_to_regress == ''){
            vars_to_regress='NONE'
        }
        
        """
        echo 'running1'
        integrate.R ${vars_to_regress} ${k_anchor} ${dims} ${ndim_sct} ${ndim_citeBgRemoved} ${ndim_cite_integrated}
        """

}

process MULTIMODAL_INTEGRATION{

    label 'process_medium'
    tag "${sample_name}"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }

    publishDir  path: "${params.outdir}/citeseq/all_data_integrated", mode: "${params.copy_mode}",
      overwrite: "true"
    
    output:
        path("*wnn.integrated.RDS"), emit: wnn_integrated_file

    input:
        path(tmp_rds_file)

    script:
    """
    echo 'running1'
    WNN_integrate_SCT_CITE.R ${tmp_rds_file}
    """

}

process VDJ_INTEGRATION{

    label 'process_medium'
    tag "${sample_name}"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }

    publishDir  path: "${params.outdir}/citeseq/all_data_integrated", mode: "${params.copy_mode}",
      overwrite: "true"
    
    output:
        path("*all_samples_integrated.vdj.RDS"), emit: all_data_integrated_vdj_rds
        path("*all_samples_integrated.BCR.RDS"), emit: all_data_integrated_BCR_rds
        path("*all_samples_integrated.TCR.RDS"), emit: all_data_integrated_TCR_rds

    input:
        path(all_cellranger_samples)
        path(wnn_integrated_file)


    script:
    """
    echo 'running1'
    add_vdj.R ${wnn_integrated_file}
    """

}


process PREPROCESS_PROCESS {
    label 'process_medium'
    tag "${sample_name}"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }

    publishDir  path: "${params.outdir}/citeseq/DSB/${sample_name}",      
        saveAs: {filename ->
        if (filename.contains("tmp_rds_files__")) {
            null
        }else {
            filename
        }
      }, mode: "${params.copy_mode}",
      overwrite: "true"

    input:
        tuple val(sample_name), path(vireo_path), path(rds_path),path(matched_donors)
        each vars_to_regress

    output:
        path("normalised__${sample_name}.withADT.RDS"), emit: tmp_rsd

    script:


        if (vars_to_regress == ''){
            vars_to_regress='NONE'
        }
        """
            process_donor_data_for_integration.R ${sample_name} ${vireo_path} ${matched_donors} ${rds_path} ${vars_to_regress}
        """
}


process DSB_PROCESS {
    label 'process_medium'
    tag "${sample_name}"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }

    publishDir  path: "${params.outdir}/citeseq/DSB/${sample_name}",      
        saveAs: {filename ->
        if (filename.contains("tmp_rds_files__")) {
            null
        }else {
            filename
        }
      }, mode: "${params.copy_mode}",
      overwrite: "true"

    input:
        tuple val(sample_name),path(cellranger_rawfile_path), path(filtered_feature_bc_matrix), path(sample_QCd_adata)
    


    output:
        path("CITE__*"), emit: citeseq_rsd
        path("tmp_rds_files__*/*/${sample_name}*.RDS"), emit: tmp_rsd
        tuple val(sample_name), path("tmp_rds_files__*/*/${sample_name}*.RDS"), emit: ch_for_norm
    script:
        """
   
            add_adt.R ${sample_name} ${cellranger_rawfile_path} ${filtered_feature_bc_matrix} ${sample_QCd_adata}
        """
}