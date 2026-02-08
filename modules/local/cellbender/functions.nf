#!/usr/bin/env nextflow

def random_hex(n) {
  Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}


if (binding.hasVariable("echo_mode") == false) {
  echo_mode = true
}

process CELLBENDER__RB__GET_INPUT_CELLS {

  label 'process_low'
  
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container "${params.yascp_container}"
    maxRetries = 1

    
  } else {
    container "${params.yascp_container_docker}"
  }

  // Calculates thresholds for input cells of CELLBENDER__REMOVE_BACKGROUND
  // ------------------------------------------------------------------------
 // use tmp directory

  publishDir  "${outdir}",
      mode: "${params.cellsnp.copy_mode}",
      overwrite: "true"

  input:
    val(outdir_prev)
    tuple(
    val(experiment_id),
    path(file_10x_barcodes),
    path(file_10x_features),
    path(file_10x_matrix),path(gex),val(droplets_include),val(expected_ncells)
    )
    val(estimate_params_umis)
    

  output:
    val(outdir, emit: outdir)
    tuple(
    val(experiment_id),
    path(file_10x_barcodes),
    path(file_10x_features),
    path(file_10x_matrix),
    path("${outfile}-expected_cells.txt"),
    path("${outfile}-total_droplets_included.txt"),
    emit: cb_input
    )
    path(
    "${outfile}-expected_cells.txt",
    emit: expected_cells
    )
    path(
    "${outfile}-total_droplets_included.txt",
    emit: total_droplets_include
    )
    path("${outfile}-cell_estimate_cutoff.tsv.gz")
    path("${outfile}-total_droplets_cutoff.tsv.gz")
    path("plots/*.png") optional true
    path("plots/*.pdf") optional true

  script:

    outdir = "${outdir_prev}/${experiment_id}"
    outdir = "${outdir}/cellbender-estimate_ncells_nemptydroplets"
    outfile = "umi_count_estimates"

    if ("${expected_ncells}"=="0"){
        cmd__expected_ncells = ""
    }else{
        cmd__expected_ncells = "--expected_ncells ${expected_ncells}"
    }

    if ("${droplets_include}"=="0"){
        cmd__droplets_include = ""
    }else{
        cmd__droplets_include = "--expected_nemptydroplets ${droplets_include}"
    }

    
    
    if ("${outdir}" == "${outdir_prev}/${experiment_id}/cellbender-params_setbyuser") {
        outdir = "${outdir_prev}/${experiment_id}/cellbender-estimated_ncells_nemptydroplets"
    }

    """
    echo ${expected_ncells}
    echo ${droplets_include}

    rm -fr plots
    mkdir txd_input
    ln --physical ${file_10x_barcodes} txd_input/barcodes.tsv.gz
    ln --physical ${file_10x_features} txd_input/features.tsv.gz
    ln --physical ${file_10x_matrix} txd_input/matrix.mtx.gz

    get_estimates_from_umi_counts.py --tenxdata_path txd_input --output_file ${outfile} --expected_nemptydroplets_umi_cutoff   ${estimate_params_umis.expected_nemptydroplets_umi_cutoff} --method_estimate_ncells   ${estimate_params_umis.method_estimate_ncells} --lower_bound_umis_estimate_ncells   ${estimate_params_umis.lower_bound_umis_estimate_ncells} --method_estimate_nemptydroplets   ${estimate_params_umis.method_estimate_nemptydroplets} --lower_bound_umis_estimate_nemptydroplets   ${estimate_params_umis.lower_bound_umis_estimate_nemptydroplets} --upper_bound_umis_estimate_nemptydroplets   ${estimate_params_umis.upper_bound_umis_estimate_nemptydroplets} --estimate_nemptydroplets_add_umifactor   ${estimate_params_umis.estimate_nemptydroplets_umi_add_factor} --estimate_nemptydroplets_subtract_dropletfactor   ${estimate_params_umis.estimate_nemptydroplets_subtract_cell_factor} --estimate_nemptydroplets_min_nemptydroplets ${estimate_params_umis.estimate_nemptydroplets_min_drop}    ${cmd__expected_ncells} ${cmd__droplets_include}

    mkdir plots
    mv *pdf plots/ 2>/dev/null || true
    mv *png plots/ 2>/dev/null || true
    """
}


process CELLBENDER__PREPROCESS_OUTPUT{
    label 'process_low'
    tag "${experiment_id}_cb"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "${params.yascp_container}"
      maxRetries = 1
      cpus = 1
    } else {
      container "${params.yascp_container_docker}"
    }
    publishDir  path: "${outdir}",
        saveAs: {filename ->
          if (filename == "barcodes.tsv.gz") {
          null
          } else if(filename.equalsIgnoreCase("features.tsv.gz")){
            null
          }
          else if(filename.equalsIgnoreCase("matrix.mtx.gz")) {
          null
          } else {
          filename.replaceAll("-", "")
          }
        },
        mode: "${params.cellsnp.copy_mode}",
        overwrite: "true" 
 
  input:

    tuple(
      val(experiment_id),
      val(fpr),
      val(outdir),
      path(filtered),
      val(cb_params)
    )

    tuple(
      val(outdir),
      path(file_10x_barcodes),
      path(file_10x_features),
      path(file_10x_matrix)
      
    )

    tuple(
      val(experiment_id),
      val(outdir),
      path(expected_cells),
      path(total_droplets_include)
    )



  output:
    tuple(val(experiment_id),path("cellbender-FPR_${params.cellbender_resolution_to_use}-filtered_10x_mtx"), emit: alternative_input)
    tuple(val(experiment_id),path("cellbender-FPR_${params.cellbender_resolution_to_use}-unfiltered_10x_mtx"), emit: alternative_input_raw)
    path("*filtered_10x_mtx/barcodes.tsv.gz", emit: tenx_barcodes)
    path("*filtered_10x_mtx/features.tsv.gz", emit: tenx_features)
    path("*filtered_10x_mtx/matrix.mtx.gz", emit: tenx_matrix)
    path("Warnings.log") optional true
    path(
      "${outfile}-filtered_10x_mtx-file_list.tsv",
       emit: results_list
    )  
    tuple(
      val(experiment_id),
      val(outdir),
      path("*_unfiltered.h5"),
      path(expected_cells),
      path(total_droplets_include),
      emit: experimentid_outdir_cellbenderunfiltered_expectedcells_totaldropletsinclude
    )
    tuple(
      val(outdir),
      path(file_10x_barcodes),
      path(file_10x_features),
      path(file_10x_matrix),
      path("*_filtered.h5"),
      emit: cb_plot_input
    )

    tuple(
      val(experiment_id),
      val(outdir),
      emit: out_paths
    )

  script:
    outfile = "cellbender"
    """
      export LD_PRELOAD=/opt/conda/envs/conda_cellbender/lib/libmkl_core.so:/opt/conda/envs/conda_cellbender/lib/libmkl_sequential.so
      for i in \$(ls *.h5); do
        echo \$i
        out_file=\$(echo \$i | sed s/"\\."/"pt"/g | sed s/"pth5"/"\\.h5"/)
        # If outfile does not exist, move i to out_file
        [ ! -f \$out_file ] && mv \$i \$out_file
        echo \$out_file
      done
      clean_cellbender_results.py --nf_outdir_tag ${outdir} --cb_outfile_tag ${outfile} --experiment_id ${experiment_id} --fpr '${fpr}' --cb_params ${cb_params}
      cp ${outfile}-filtered_10x_mtx-file_list.tsv ${outfile}-filtered_10x_mtx-file_list.tsv || echo 'same file'
    """

}

process CELLBENDER__REMOVE_BACKGROUND {
  // Remove ambient RNA
  // ------------------------------------------------------------------------
  //cache false    // cache results from run
  //maxForks 2   // hard to control memory usage. limit to 3 concurrent
// cb_plot_input,out_paths,results_list,experimentid_outdir_cellbenderunfiltered_expectedcells_totaldropletsinclude

  tag "${experiment_id}_cb"
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    if (params.cellbender_v == '0.3.2'){
      container "${params.nf_cellbender_container_032}"
    }else{
      container "${params.nf_cellbender_container}"
    }
    
    maxRetries = 1
	  cpus = 1

  } else {
    if (params.cellbender_v == '0.3.2'){
      container "${params.nf_cellbender_container_032_docker}"
    }else{
      container "${params.nf_cellbender_container_docker}"
    }
    
  }
    
  // set LD_PRELOAD to fix mkl / anaconda conflict
  // cf. https://stackoverflow.com/questions/36659453/intel-mkl-fatal-error-cannot-load-libmkl-avx2-so-or-libmkl-def-so
  // containerOptions '--env LD_PRELOAD=/opt/conda/envs/conda_cellbender/lib/libmkl_core.so:/opt/conda/envs/conda_cellbender/lib/libmkl_sequential.so'

  //     // use GPU
  if (params.utilise_gpu){
    label 'gpu'

	// only one label here, otherwise bsub -R -M will be doubled..
	// label 'process_high_memory'
    gpu_text_info = '--cuda'
  }else{
    label 'process_medium'
  }

  // scratch false    // use tmp directory
  // echo false   // echo output from script

  publishDir  path: "${outdir}",
      saveAs: {filename ->
        if (filename == "barcodes.tsv.gz") {
        null
        } else if(filename.equalsIgnoreCase("features.tsv.gz")) {
        null
        } else if(filename.equalsIgnoreCase("matrix.mtx.gz")) {
        null
        } else {
        filename.replaceAll("-", "")
        }
      },
      mode: "${params.cellsnp.copy_mode}",
      overwrite: "true"

  input:
    val(outdir_prev)
    tuple(
      val(experiment_id),
      path(file_10x_barcodes),
      path(file_10x_features),
      path(file_10x_matrix),
      path(expected_cells),
      path(total_droplets_include),
      val(low_count_threshold),
      val(epochs),
      val(learning_rate),
      val(zdims),
      val(zlayers),
    )
    // each epochs
    // each learning_rate
    // each zdims
    // each zlayers
    // each low_count_threshold
    each fpr

  output:
    val(outdir, emit: outdir)
    tuple(
      val(experiment_id),
      val(fpr),
      val(outdir),
      path("*.h5"),val(cb_params),
      emit: cleanup_input
    )
    tuple(val(experiment_id),path("cellbender_FPR_${params.cellbender_resolution_to_use.replaceAll('pt', '.')}_filtered.h5"),emit: cb_to_use_downstream)

    tuple(
      val(outdir),
      val(experiment_id),
      val(epochs),
      val(learning_rate),
      val(zdims),
      val(zlayers),
      val(low_count_threshold),
      val(fpr),
      emit: cb_results_details
    )

    tuple(
      val(outdir),
      path(file_10x_barcodes),
      path(file_10x_features),
      path(file_10x_matrix),
      emit: cb_plot_input
    )

    path("${outfile}_cell_barcodes.csv", emit: barcodes)
    
    path("*.h5", emit: filtered_h5s)
    path(file_10x_barcodes, emit: raw_tenx_barcodes)
    path(file_10x_features, emit: raw_tenx_features)
    path(file_10x_matrix, emit: raw_tenx_matrix)
    path("plots/*.png") optional true
    path("plots/*.pdf") optional true

    tuple(
      val(experiment_id),
      val(outdir),
      
      path(expected_cells),
      path(total_droplets_include),
      emit: experimentid_outdir_cellbenderunfiltered_expectedcells_totaldropletsinclude
    )

    tuple(
      val(experiment_id),
      val(outdir),
      emit: out_paths
    )
    path "versions.yml", emit: versions

  script:

    if (params.cellbender_v == '0.3.1'){
        option1='--checkpoint-mins 100'
        if (params.utilise_gpu){
            gpu_text_info = '--cuda'
        }else{
            gpu_text_info = "--cpu-threads ${task.cpus}"
        }
    }else{
        option1=''
        if (params.utilise_gpu){
            gpu_text_info = '--cuda'
        }else{
            gpu_text_info = ""
        }
    }

    outdir = "${outdir_prev}/${experiment_id}"
    lr_string = "${learning_rate}".replaceAll("\\.", "pt")
    lr_string = "${lr_string}".replaceAll("-", "neg")
    fpr_string = "${fpr}".replaceAll("\\.", "pt").replaceAll(" ", "_")
    cb_params = "cellbender_params"
    cb_params = "${cb_params}-epochs_${epochs}"
    cb_params = "${cb_params}__learnrt_${lr_string}"
    cb_params = "${cb_params}__zdim_${zdims}"
    cb_params = "${cb_params}__zlayer_${zlayers}"
    cb_params = "${cb_params}__lowcount_${low_count_threshold}"
    outdir = "${outdir}/${cb_params}".replaceAll("cellbender_params","cellbender")
    outfile = "cellbender"
  
      """
    echo ${experiment_id}

    export LD_PRELOAD=/opt/conda/envs/conda_cellbender/lib/libmkl_core.so:/opt/conda/envs/conda_cellbender/lib/libmkl_sequential.so

    rm -fr plots
    mkdir -p txd_input

    # Check if input barcodes are already unique 16bp identifiers.
    # If not, this script will generate deterministic 16-character encoded barcodes (CB00000000000000, etc.),
    # write them to txd_input/barcodes_encoded.tsv.gz, and create a barcode mapping file
    # at txd_input/barcode_mapping.tsv for downstream remapping of CellBender outputs.

    encode_barcodes_if_needed.py

    # Link features and matrix
    ln --physical ${file_10x_features} txd_input/features.tsv.gz
    ln --physical ${file_10x_matrix} txd_input/matrix.mtx.gz

    # Conditionally link encoded or original barcodes
    if [[ -f txd_input/barcodes_encoded.tsv.gz ]]; then
        echo "Using encoded barcodes"
        ln --symbolic barcodes_encoded.tsv.gz txd_input/barcodes.tsv.gz
    else
        echo "Using original barcodes"
        ln --physical ${file_10x_barcodes} txd_input/barcodes.tsv.gz
    fi

    export TMPDIR=\$PWD
    cellbender remove-background \\
      --input txd_input \\
      ${gpu_text_info} ${option1} \\
      --output ${outfile} \\
      --expected-cells \$(cat ${expected_cells}) \\
      --total-droplets-included \$(cat ${total_droplets_include}) \\
      --model full \\
      --z-dim ${zdims} \\
      --z-layers ${zlayers} \\
      --low-count-threshold ${low_count_threshold} \\
      --epochs ${epochs} \\
      --learning-rate ${learning_rate} \\
      --fpr ${fpr}

    # Rename output if needed
    [ -f ${outfile} ] && mv ${outfile} ${outfile}.h5

    # Remap h5 barcodes if mapping file exists
    if [[ -f txd_input/barcode_mapping.tsv ]]; then
        echo "Remapping .h5 barcodes"
        resave_h5_with_mapped_barcodes.py
    fi

    mkdir -p plots
    mv *pdf plots/ 2>/dev/null || true
    mv *png plots/ 2>/dev/null || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellbender: \$(pip show cellbender | grep '^Version:' | cut -d':' -f2)
    END_VERSIONS
    """
}

process CAPTURE_CELLBENDER_FILES{

  publishDir  path: "${outdir}/cellbender",
  saveAs: {filename ->
        if (filename.contains("captured")) {
          null
        } else if (filename.contains(".h5")){
          null
        }
        else {
          filename.replaceAll("tmp1234/cellbender/", "")
        }
      },
        mode: "${params.copy_mode}",
    overwrite: "true"
  label 'process_tiny'



  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    container "${params.yascp_container}"
  } else {
    container "${params.yascp_container_docker}"
  }

  input:
    path(cellbender_location) //cant use path here - needs a full path
    val(outdir)
    path(input_file)
  output:
    path("tmp1234/cellbender/*") optional true
    path("captured/unfiltered/*/*FPR_${params.cellbender_resolution_to_use}*"),emit:alt_input_unfiltered optional true
    path("captured/filtered/*/*FPR_${params.cellbender_resolution_to_use}*"),emit:alt_input_filtered optional true

  script:
  """
  
  
    mkdir tmp1234
    mkdir tmp1234/cellbender
    cd tmp1234/cellbender
    for d in ../../${cellbender_location}/* ; do
        echo \$d
        basename1=\$(basename \$d)
        if grep -q \$basename1 ../../${input_file}
        then
          ln -s \$d .
        else
            f='no'
        fi        
    done
    cd ../..
    capture_res_files_cb.py -res ${params.cellbender_resolution_to_use}
  """    

}

process CELLBENDER__REMOVE_BACKGROUND__QC_PLOTS {

  label 'process_low'
  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"

    } else {
        container "${params.yascp_container_docker}"
  }

  // Second set of QC plots from cellbdender.
  // This task compare the cellbender output with both the cellranger filtered and cellragner raw outputs
  // ------------------------------------------------------------------------
  tag { "$experiment_id" }
  publishDir "${outdir2}/$experiment_id/compare_cellranger/", pattern: "fpr_${fpr}/${experiment_id}/*.png", 
    saveAs: {filename ->
    filename.replaceAll("fpr_${fpr}/${experiment_id}/", "fpr_${fpr}/")
  },
    mode: "${params.cellsnp.copy_mode}",
    overwrite: "true"
  
  //cache false    // cache results from run
  //maxForks 2   // hard to control memory usage. limit to 3 concurrent
  // scratch false    // use tmp directory
  // echo echo_mode   // echo output from script
  input:
    tuple val(experiment_id), val(outdir), path(cellbender_unfiltered_h5s), path(expectedcells), path(totaldropletsinclude), path(raw_cellranger_mtx), path(filtered_cellranger_mtx), val(fpr)
    val(outdir2)
  output:
    val(outdir, emit: outdir)
    path("fpr_${fpr}/${experiment_id}/*.png"), emit: plots_png 

  script:
  """
  fprid=\$(echo $fpr | sed s'/\\./pt/'g)
  cellbender_unfiltered_h5=cellbender_FPR_\${fprid}_unfiltered.h5

  n_expected_cells=\$(cat $expectedcells)
  n_total_droplets_included=\$(cat $totaldropletsinclude)

  plot_cellranger_vs_cellbender.py \\
    --samplename \"${experiment_id}\" \\
    --raw_cellranger_mtx \"${raw_cellranger_mtx}\" \\
    --filtered_cellranger_mtx \"${filtered_cellranger_mtx}\" \\
    --cellbender_unfiltered_h5 \"\$cellbender_unfiltered_h5\" \\
    --fpr \"${fpr}\" \\
    --n_expected_cells \"\${n_expected_cells}\" \\
    --n_total_droplets_included \"\${n_total_droplets_included}\" \\
    --out_dir \$PWD
  """
}

