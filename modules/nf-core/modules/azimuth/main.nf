process AZIMUTH{
    tag "${samplename}"    
    label 'process_medium'
   
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://yascp.cog.sanger.ac.uk/public/singularity_images/wtsihgi_nf_scrna_qc_azimuth_d54db9b-2021-12-13-8dd0b7fce918.sif"
        //// container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/seurat_azimuth_pbmc_1.0.img"
        
    } else {
        container "wtsihgi/nf_scrna_qc_azimuth:d54db9b"
    }

    publishDir  path: "${params.outdir}/celltype/azimuth/${refset.name}",
            saveAs: {filename -> "${outfil_prfx}_" + filename},
            mode: "${params.copy_mode}",
            overwrite: "true"
    stageInMode 'copy'
    // stageInMode 'copy' is needed because SeuratDisk:::Convert()
    // generates the output file apparently from the absolute path of input file.
    // Symbolic links have the output file written to the link target directory
    // where it cannot be found by the azimuth.R script.

    input:
        val outdir_prev
        tuple val(samplename),path(file_h5ad_batch)
        each refset
        // path(mapping_file)
    output:
        // path(celltype_table, emit:predicted_celltypes)
        tuple(val(outfil_prfx), val(refset.refset), path("*predicted_*.tsv"),emit:celltype_tables_all) 
        path("*predicted_*.tsv"), emit:predicted_celltype_labels
        path "*ncells_by_type_barplot.pdf"
        path "*query_umap.pdf"
        path "*prediction_score_umap.pdf"
        path "*prediction_score_vln.pdf"
        path "*mapping_score_umap.pdf"
        path "*mapping_score_vln.pdf"
        // path("${outfil_prfx}_query.rds"), emit: query_rds

        
    script:
    
    
    // output file prefix: strip random hex number form beginning of file name
    outfil_prfx = "${file_h5ad_batch}".minus(".h5ad")
    //outfil_prfx = "${file_h5ad_batch}".minus(".h5ad")
    
    """ 
        azimuth.R ./${file_h5ad_batch} ${refset.refset} ${refset.annotation_labels} ${samplename}

    """
}

process REMAP_AZIMUTH{
    // This process remaps Azimuth L2 to L1 and L0
    tag "${samplename}"    
    label 'process_low'
   
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://yascp.cog.sanger.ac.uk/public/singularity_images/wtsihgi_nf_scrna_qc_6bb6af5-2021-12-23-3270149cf265.sif"
        //// container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/seurat_azimuth_pbmc_1.0.img"
    } else {
        container "wtsihgi/nf_scrna_qc_azimuth:d54db9b"
    }

    publishDir  path: "${params.outdir}/celltype/azimuth/",
            mode: "${params.copy_mode}",
            overwrite: "true"
    stageInMode 'copy'  

    input:
        tuple val(outfil_prfx),val(name) path(azimuth_file)
        path(mapping_file)
    when:
        name=='PBMC'
    output:
        path('azimuth/*', emit:predicted_celltype_labels)

    script:
        celltype_table = "remapped__${azimuth_file}"
        """
            remap_azimuth_l2.py -of remapped__predicted_celltype_l2.tsv -m ${mapping_file} -az predicted_celltype_l2.tsv
            mkdir azimuth
            
            mv predicted_celltype_l1.tsv  azimuth/${outfil_prfx}__predicted_celltype_l1.tsv
            mv predicted_celltype_l3.tsv azimuth/${outfil_prfx}__predicted_celltype_l3.tsv
            mv remapped__predicted_celltype_l2.tsv azimuth/${outfil_prfx}__predicted_celltype_l2.tsv
        """

}
