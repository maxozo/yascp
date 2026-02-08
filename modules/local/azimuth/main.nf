process AZIMUTH{
    tag "${samplename}"    
    label 'process_medium'
   
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"        
    } else {
        container "${params.yascp_container_docker}"
    }

    publishDir  path: "${params.outdir}/celltype_assignment/azimuth/${refset.name}",
            saveAs: {filename -> "${outfil_prfx}_" + filename},
            mode: "${params.copy_mode}",
            overwrite: "true"
    stageInMode 'copy'
    // stageInMode 'copy' is needed because SeuratDisk:::Convert()
    // generates the output file apparently from the absolute path of input file.
    // Symbolic links have the output file written to the link target directory
    // where it cannot be found by the azimuth.R script.

    input:
        tuple val(samplename),path(file_h5ad_batch)
        each path(mapping_file)
        each refset
    output:
        tuple(val(outfil_prfx), val(refset.refset), path("*predicted_*.tsv"),emit:celltype_tables_all) 
        path("*predicted_*.tsv"), emit:predicted_celltype_labels
        path "*ncells_by_type_barplot.pdf"
        path "*query_umap.pdf"
        path "*prediction_score_umap.pdf"
        path "*prediction_score_vln.pdf"
        path "*mapping_score_umap.pdf"
        path "*mapping_score_vln.pdf"
        path "versions.yml", emit: versions

    script:

        // output file prefix: strip random hex number form beginning of file name
        outfil_prfx = "${file_h5ad_batch}".minus(".h5ad")
        //outfil_prfx = "${file_h5ad_batch}".minus(".h5ad")
        if (refset.refset =='PBMC' && params.mapping_file!='' && params.remap_celltypes){
            com="remap_azimuth_l2.R --out_file ${samplename}___predicted_celltype_l2.tsv --mapping ${mapping_file} --az_file ${samplename}___predicted_celltype_l2.tsv"
        }else{
            com=""
        }

        if (params.atac){
            atac = '--atac'
        }else{
            atac = ""
        }    

    """ 
        azimuth.R ./${file_h5ad_batch} ${refset.refset} ${refset.annotation_labels} ${samplename}
        ${com}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            r-base: \$(R --version | sed -n '1p' | sed 's/R version //; s/ (.*//')
            r library Azimuth: \$(Rscript -e "cat(as.character(packageVersion('Azimuth')))")
            r library ggplot2: \$(Rscript -e "cat(as.character(packageVersion('ggplot2')))")
            r library hdf5r: \$(Rscript -e "cat(as.character(packageVersion('hdf5r')))")
            r library Matrix: \$(Rscript -e "cat(as.character(packageVersion('Matrix')))")
            r library Seurat: \$(Rscript -e "cat(as.character(packageVersion('Seurat')))")
            r library SeuratDisk: \$(Rscript -e "cat(as.character(packageVersion('SeuratDisk')))")
            r library tools: \$(Rscript -e "cat(as.character(packageVersion('tools')))")
        END_VERSIONS

    """
}


process AZIMUTH_ATAC{
    tag "${samplename}"    
    label 'process_medium'
   
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"        
    } else {
        container "${params.yascp_container_docker}"
    }

    publishDir  path: "${params.outdir}/celltype_assignment/azimuth/${refset.name}",
            saveAs: {filename -> "${outfil_prfx}_" + filename},
            mode: "${params.copy_mode}",
            overwrite: "true"
    input:
        tuple val(samplename),path(file_h5ad_batch)
        each path(mapping_file)
        each refset

    script: 

    """ 
        azimuth_atac.R ./${file_h5ad_batch} ${refset.refset} ${refset.annotation_labels} ${samplename}
    """
}
