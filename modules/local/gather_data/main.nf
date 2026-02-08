process GATHER_DATA{

    publishDir  path: "${params.outdir}/handover",
                mode: "${params.copy_mode}",
                overwrite: "true"
    label 'process_medium'

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }

    input:
      path(outdir_prev)
      val(dummy_val)
      path(input_data_table)

    output:
      path("${subdir}", emit:outfiles_dataset) optional true
      path("${subdir}_summary", emit:outfiles_dataset2) optional true
      path("Donor_Quantification/*/*.tsv", emit: barcodes_files) optional true
      val(outdir, emit: outdir_dataset)

    script:
      outdir = "${outdir_prev}/handover"
      subdir = "Donor_Quantification"
      if ("${params.input}" == 'cellranger'){
        cellbender_input='cellranger'
      }else{
        cellbender_input='cellbender'
      }

      if ("${params.extra_sample_metadata}" != ''){
        extra_meta = "--extra_meta=${params.extra_sample_metadata}"
      }else{
        extra_meta = ""
      }

      """
        echo ${dummy_val}
        gather_minimal_dataset.py \
          --output-dir=${subdir} \
          --results_dir=${outdir_prev} \
          --input_table=${input_data_table} \
          --cellbender=${cellbender_input} \
          --resolution=${params.cellbender_resolution_to_use} \
          --write_h5=False \
          --experiment_name=${params.RUN} ${extra_meta}
      """
}
