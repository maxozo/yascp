process DOUBLET_FINDER {

    tag "${experiment_id}"
    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }
    
    publishDir  path: "${params.outdir}/doublet_detection/DoubletFinder",
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        tuple(
            val(experiment_id),
            path(h5ad)
        )
        val(expected_multiplet_rate)

    output:
        path("plots/*.pdf") optional true
        path("plots/*.png") optional true
        tuple val(experiment_id), path("${experiment_id}__DoubletFinder_doublets_singlets.tsv"), emit: result optional true
        path "versions.yml", emit: versions

    script:
        
        outdir = "${params.outdir}/doublet_detection/"
        outdir = "${outdir}doubletFinder"
        outfile = "${experiment_id}"

        """
            DoubletFinder.R -s ${h5ad} -c FALSE -o DoubletFinder_${experiment_id} --expected_multiplet_rate ${expected_multiplet_rate}
            if [ -e "DoubletFinder_${experiment_id}/DoubletFinder_doublets_singlets.tsv" ]; then
                ln -s "DoubletFinder_${experiment_id}/DoubletFinder_doublets_singlets.tsv" "${experiment_id}__DoubletFinder_doublets_singlets.tsv"
            else
                echo "File DoubletFinder_${experiment_id}/DoubletFinder_doublets_singlets.tsv does not exist. Skipping symlink creation."
            fi

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            r-base: \$(R --version | sed -n '1p' | sed 's/R version //; s/ (.*//')
            r library argparse: \$(Rscript -e "cat(as.character(packageVersion('argparse')))")
            r library DoubletFinder: \$(Rscript -e "cat(as.character(packageVersion('DoubletFinder')))")
            r library dplyr: \$(Rscript -e "cat(as.character(packageVersion('dplyr')))")
            r library future.apply: \$(Rscript -e "cat(as.character(packageVersion('future.apply')))")
            r library ggplot2: \$(Rscript -e "cat(as.character(packageVersion('ggplot2')))")
            r library Seurat: \$(Rscript -e "cat(as.character(packageVersion('Seurat')))")
            r library SeuratDisk: \$(Rscript -e "cat(as.character(packageVersion('SeuratDisk')))")
            r library tidyr: \$(Rscript -e "cat(as.character(packageVersion('tidyr')))")
            r library tidyverse: \$(Rscript -e "cat(as.character(packageVersion('tidyverse')))")
            r library viridis: \$(Rscript -e "cat(as.character(packageVersion('viridis')))")
        END_VERSIONS
        """
}
