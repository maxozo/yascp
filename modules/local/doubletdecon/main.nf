process DOUBLET_DECON{

    tag "${experiment_id}"
    // label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }
    
    publishDir  path: "${params.outdir}/doublet_detection/DoubletDecon",
                mode: "${params.copy_mode}",
                overwrite: "true"
    memory { 
            sizeInGB = h5ad.size() / 1e9 * 5 * task.attempt
            return (sizeInGB ).toString() + 'GB' 
        }
        
    input:
        tuple(
            val(experiment_id),
            path(h5ad)
        )

    output:
        path("plots/*.pdf") optional true
        path("plots/*.png") optional true
        tuple val(experiment_id), path("${experiment_id}__DoubletDecon_doublets_singlets.tsv"), emit: result  optional true
        path "versions.yml", emit: versions

    script:
        
        outdir = "${params.outdir}/doublet_detection/"
        outdir = "${outdir}DoubletDecon"
        outfile = "${experiment_id}"

        """
            DoubletDecon.R -s ${h5ad} -o DoubletDecon_${experiment_id}
            
            if [ -e "DoubletDecon_${experiment_id}/DoubletDecon_doublets_singlets.tsv" ]; then
                ln -s DoubletDecon_${experiment_id}/DoubletDecon_doublets_singlets.tsv ${experiment_id}__DoubletDecon_doublets_singlets.tsv
            else
                echo "File DoubletDecon_${experiment_id}/DoubletDecon_doublets_singlets.tsv does not exist. Skipping symlink creation."
            fi

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                r-base: \$(R --version | sed -n '1p' | sed 's/R version //; s/ (.*//')
                r library argparse: \$(Rscript -e "cat(as.character(packageVersion('argparse')))")
                r library data.table: \$(Rscript -e "cat(as.character(packageVersion('data.table')))")
                r library DoubletDecon: \$(Rscript -e "cat(as.character(packageVersion('DoubletDecon')))")
                r library ggplot2: \$(Rscript -e "cat(as.character(packageVersion('ggplot2')))")
                r library future.apply: \$(Rscript -e "cat(as.character(packageVersion('future.apply')))")
                r library Seurat: \$(Rscript -e "cat(as.character(packageVersion('Seurat')))")
                r library SeuratDisk: \$(Rscript -e "cat(as.character(packageVersion('SeuratDisk')))")
                r library tidyverse: \$(Rscript -e "cat(as.character(packageVersion('tidyverse')))")
                r library viridis: \$(Rscript -e "cat(as.character(packageVersion('viridis')))")
            END_VERSIONS
        """
}
