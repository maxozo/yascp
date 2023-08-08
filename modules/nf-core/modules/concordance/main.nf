process CONCORDANCE_CALCLULATIONS {

    tag "${pool_id}"
    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://yascp.cog.sanger.ac.uk/public/singularity_images/mercury_scrna_deconvolution_62bd56a-2021-12-15-4d1ec9312485.sif"
        //// container "https://yascp.cog.sanger.ac.uk/public/singularity_images/mercury_scrna_deconvolution_latest.img"
    } else {
        container "mercury/scrna_deconvolution:62bd56a"
    }

    publishDir  path: "${params.outdir}/concordances/${pool_id}",
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        tuple(val(pool_id), 
        path(vcf_gt_match), 
        path(vcf_gt_match_csi),
        path(vcf_exp), 
        path(vcf_exp_csi),
        path(cell_vcf),
        path(donor_table),path(cell_assignments),path(set2_informative_sites), path(set1_uninformative_sites),path(variants_description))

    output:
        tuple val(pool_id), path("cell_concordance_table_noA2G.tsv"), emit: concordances_noAG
        tuple val(pool_id), path('discordant_sites_in_other_donors_noA2G.txt'), emit: read_concordances_noAG
        tuple val(pool_id), path("cell_concordance_table.tsv"), emit: concordances
        tuple val(pool_id), path('discordant_sites_in_other_donors.txt'), emit: read_concordances

    script:

        """
            echo ${pool_id}
            bcftools view ${vcf_exp} -R ${cell_vcf} -Oz -o sub_${pool_id}_Expected.vcf.gz
            if [ "${vcf_gt_match}" != "fake_file.fq" ]; then
                bcftools view ${vcf_gt_match} -R ${cell_vcf} -Oz -o sub_${pool_id}_GT_Matched.vcf.gz
            fi
            
            concordance_calculations_donor_exclusive_read_level_noA2G.py --cpus $task.cpus --cell_vcf ${cell_vcf} --donor_assignments ${donor_table} --gt_match_vcf sub_${pool_id}_GT_Matched.vcf.gz --expected_vcf sub_${pool_id}_Expected.vcf.gz --cell_assignments ${cell_assignments} --informative_sites ${set2_informative_sites} --uninformative_sites ${set1_uninformative_sites}
            find_discordant_sites_in_other_donors_noA2G.py --cpus $task.cpus --cell_vcf ${cell_vcf} --donor_assignments ${donor_table} --gt_match_vcf sub_${pool_id}_GT_Matched.vcf.gz --expected_vcf sub_${pool_id}_Expected.vcf.gz --cell_assignments ${cell_assignments} --outfile discordant_sites_in_other_donors_noA2G.txt --debug
            concordance_calculations_donor_exclusive_read_level.py --cpus $task.cpus --cell_vcf ${cell_vcf} --donor_assignments ${donor_table} --gt_match_vcf sub_${pool_id}_GT_Matched.vcf.gz --expected_vcf sub_${pool_id}_Expected.vcf.gz --cell_assignments ${cell_assignments} --informative_sites ${set2_informative_sites} --uninformative_sites ${set1_uninformative_sites}
            find_discordant_sites_in_other_donors_find_best_donor.py --cpus $task.cpus --cell_vcf ${cell_vcf} --donor_assignments ${donor_table} --gt_match_vcf sub_${pool_id}_GT_Matched.vcf.gz --expected_vcf sub_${pool_id}_Expected.vcf.gz --cell_assignments ${cell_assignments} --outfile discordant_sites_in_other_donors.txt --debug
        """
}


process COMBINE_FILES{

    tag "${pool_id}"
    label 'process_low'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://yascp.cog.sanger.ac.uk/public/singularity_images/mercury_scrna_deconvolution_62bd56a-2021-12-15-4d1ec9312485.sif"
        //// container "https://yascp.cog.sanger.ac.uk/public/singularity_images/mercury_scrna_deconvolution_latest.img"
    } else {
        container "mercury/scrna_deconvolution:62bd56a"
    }

    publishDir  path: "${params.outdir}/concordances/${pool_id}",
                mode: "${params.copy_mode}",
                overwrite: "true"


    input:
        tuple(val(pool_id), 
        path(subsampling_table), 
        path(concordance_table))

    output:
        path("*.png")
        tuple val(pool_id), path("*joined_df_for_plots.tsv"), emit: joined_df_for_plots
        path("*joined_df_for_plots.tsv"), emit: file_joined_df_for_plots

    script:

        """
           combine_concordance.py -cc ${concordance_table} -sq ${subsampling_table} -name ${pool_id} --run ${params.RUN}
        """

}


process PLOT_CONCORDANCES_ALL{

    tag "${pool_id}"
    label 'process_low'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://yascp.cog.sanger.ac.uk/public/singularity_images/mercury_scrna_deconvolution_62bd56a-2021-12-15-4d1ec9312485.sif"
        //// container "https://yascp.cog.sanger.ac.uk/public/singularity_images/mercury_scrna_deconvolution_latest.img"
    } else {
        container "mercury/scrna_deconvolution:62bd56a"
    }

    publishDir  path: "${params.outdir}/concordances",
                mode: "${params.copy_mode}",
                overwrite: "true"

    output:
        path("*.png")

    input:
        path(input_file_all)

    script:

        """
            combined_concordance_plots.py -cc ${input_file_all} -name all
        """


}
