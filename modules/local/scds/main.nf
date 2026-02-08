process SCDS {

    tag "${experiment_id}"
    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }
    
    publishDir  path: "${params.outdir}/doublet_detection/SCDS",
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        tuple(
            val(experiment_id),
            path(file_10x)
        )
    output:
        path("plots/*.pdf") optional true
        path("plots/*.png") optional true
        tuple val(experiment_id), path("${experiment_id}__scds_doublets_singlets.tsv"), emit: result
        path "versions.yml", emit: versions

    script:
        
        outdir = "${params.outdir}/"
        outdir = "${outdir}SCDS"
        outfile = "${experiment_id}"

        """
            mkdir scds_${experiment_id}
            scds.R -t ${file_10x} -o scds_${experiment_id}
            cat scds_${experiment_id}/scds_doublets_singlets.tsv | awk -F' ' '{print \$1"\\t"\$3"\\t"\$2}' > ${experiment_id}__scds_doublets_singlets.tsv

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                r-base: \$(R --version | sed -n '1p' | sed 's/R version //; s/ (.*//')
                r library argparse: \$(Rscript -e "cat(as.character(packageVersion('argparse')))")
                r library dplyr: \$(Rscript -e "cat(as.character(packageVersion('dplyr')))")
                r library scds: \$(Rscript -e "cat(as.character(packageVersion('scds')))")
                r library Seurat: \$(Rscript -e "cat(as.character(packageVersion('Seurat')))")
                r library SingleCellExperiment: \$(Rscript -e "cat(as.character(packageVersion('SingleCellExperiment')))")
                r library tidyr: \$(Rscript -e "cat(as.character(packageVersion('tidyr')))")
                r library tidyverse: \$(Rscript -e "cat(as.character(packageVersion('tidyverse')))")
            END_VERSIONS
        """
}