process DOUBLET_DECON{

    tag "${experiment_id}"
    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://yascp.cog.sanger.ac.uk/public/singularity_images/azimuth_dsb_6_03_2024.sif"
    } else {
        container "mercury/azimuth_dsb:6_03_2024"
    }
    
    publishDir  path: "${params.outdir}/doublets/multiplet.method=DoubletDecon",
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        tuple(
            val(experiment_id),
            path(h5ad)
        )

    output:
        path("plots/*.pdf") optional true
        path("plots/*.png") optional true
        tuple val(experiment_id), path("${experiment_id}__DoubletDecon_doublets_singlets.tsv"), emit: result
        // path("${experiment_id}__DoubletDetection_results.txt"), emit: doubletDetection_results

    script:
        
        outdir = "${params.outdir}/doublets/multiplet"
        outdir = "${outdir}.method=DoubletDecon"
        outfile = "${experiment_id}"

        """
            DoubletDecon.R -s ${h5ad} -o DoubletDecon_${experiment_id}
            ln -s DoubletDecon_${experiment_id}/DoubletDecon_doublets_singlets.tsv ${experiment_id}__DoubletDecon_doublets_singlets.tsv
        """
}
