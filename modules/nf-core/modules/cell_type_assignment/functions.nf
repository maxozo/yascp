process CELLTYPE_FILE_MERGE{
    tag "${samplename}"    
    label 'process_high'
    publishDir  path: "${params.outdir}/celltype/",
            saveAs: {filename -> filename},
            mode: "${params.copy_mode}",
            overwrite: "true"  
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://yascp.cog.sanger.ac.uk/public/singularity_images/wtsihgi_nf_scrna_qc_6bb6af5-2021-12-23-3270149cf265.sif"
        // container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_qc_cluster_2.4.img"
    } else {
        container "wtsihgi/nf_scrna_qc:6bb6af5"
    }
    output:
        path('adata.h5ad', emit:file__anndata_merged2)
        path("All_Celltype_Assignments.tsv",emit:celltype_assignments)
        path "tranche_celltype_report.tsv"
        path "donor_celltype_report.tsv"

    input:
        path(azimuth_files)
        path(celltypist_paths)
        path(all_other_paths)
        path(file__anndata_input)
    script:
        def merged_files_outpath = "${params.outdir}/celltype/merged_files/"
        file(merged_files_outpath).mkdirs()
        def azimuth_files_path = "${merged_files_outpath}/azimuth_files.tsv"
        def celltypist_files_path = "${merged_files_outpath}/celltypist_files.tsv"
        def all_other_files_path = "${merged_files_outpath}/other_files.tsv"
        def adatas_path = "${merged_files_outpath}/adatas.tsv"

        new File(azimuth_files_path).text = azimuth_files.join("\n")
        new File(celltypist_files_path).text = celltypist_paths.join("\n")

        if ("${all_other_paths}" != 'fake_file.fq') {
            new File(all_other_files_path).text = all_other_paths.join("\n")
            other_paths = "--all_other_paths ${all_other_files_path}"
        } else {
            other_paths = ""
        }

        new File(adatas_path).text = file__anndata_input.join("\n")

        """
        generate_combined_celltype_anotation_file.py --all_azimuth_files ${azimuth_files_path} --all_celltypist_files ${celltypist_files_path} ${other_paths} --adata '${adatas_path}'
        """

}
