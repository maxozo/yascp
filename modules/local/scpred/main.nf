process SCPRED{
    tag "${samplename}"    
    label 'process_medium'
   
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }

    publishDir  path: "${params.outdir}/celltype_assignment/scpred/",
            saveAs: {filename -> "${outfil_prfx}_" + filename},
            mode: "${params.copy_mode}",
            overwrite: "true"
    stageInMode 'copy'

    input:
        tuple val(exp_id), path(file_h5ad_batch)
        path(reference)
    when:
        params.celltype_assignment.run_scpred
    output:
        path("${outfil_prfx}___scpred_prediction.tsv"), emit:predicted_celltype_labels
        path "versions.yml", emit: versions

    script:

        outfil_prfx = "${file_h5ad_batch}".minus(".h5ad")
        celltype_table = "${outfil_prfx}_predicted_celltype_l2.tsv.gz"
        """
            scpred_map_hiers.R --file ./${file_h5ad_batch} --reference ${reference}
            ln -s scpred_prediction.tsv ${outfil_prfx}___scpred_prediction.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            r-base: \$(R --version | sed -n '1p' | sed 's/R version //; s/ (.*//')
            r library future.apply: \$(Rscript -e "cat(as.character(packageVersion('future.apply')))")
            r library HierscPred: \$(Rscript -e "cat(as.character(packageVersion('HierscPred')))")
            r library optparse: \$(Rscript -e "cat(as.character(packageVersion('optparse')))")
            r library progressr: \$(Rscript -e "cat(as.character(packageVersion('progressr')))")
            r library Seurat: \$(Rscript -e "cat(as.character(packageVersion('Seurat')))")
            r library SeuratDisk: \$(Rscript -e "cat(as.character(packageVersion('SeuratDisk')))")
        END_VERSIONS
        """
}
