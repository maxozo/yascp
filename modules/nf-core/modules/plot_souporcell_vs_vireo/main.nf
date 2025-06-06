
process SOUPORCELL_VS_VIREO {

    tag "${samplename}"
    publishDir "${params.outdir}/deconvolution/souporcell_vs_vireo/",
	mode: "${params.plot_souporcell_vs_vireo.copy_mode}",
	pattern: "${samplename}_souporcell_vs_vireo.pdf",
	overwrite: true
    
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }


    when: 
    params.plot_souporcell_vs_vireo.run

    input: 
    tuple val(samplename), path(donor_ids_tsv), path(clusters_tsv)
    
    output:
    tuple val(samplename), file("${samplename}_souporcell_vs_vireo.pdf"), emit: plot_pdf

    script:
    """
        umask 2 # make files group_writable
        Rscript --vanilla ${projectDir}/../bin/plot_souporcell_vs_vireo.R ${samplename} ${donor_ids_tsv} ${clusters_tsv}
    """
}

