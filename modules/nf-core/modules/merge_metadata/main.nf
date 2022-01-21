process prep_collectmetadata{
    input:
        tuple val(experiment_id), path(metadata_path)
    output:
        path("${experiment_id}---metadata.csv", emit: metadata)
    script:
        """
            ln --physical ${metadata_path} ${experiment_id}---metadata.csv
        """
}

process merge_metadata{

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
         container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_qc_cluster_2.4.img"
     } else {
         container "quay.io/biocontainers/multiqc:1.10.1--py_0"
     }

     input:
        file(file_metadata)
    output:
        path("full_metadata.tsv", emit: metadata)
    script:
        files__metadata = file_metadata.join(',')
        """
            python ${projectDir}/bin/combine_metadata.py -d ${files__metadata}
        """
}
