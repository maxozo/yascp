def random_hex(n) {
    Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}

process PLOT_FILTERED_CELLS {
    // Takes annData object, plots filtered cells
    // ------------------------------------------------------------------------
    //cache false        // cache results from run
    label 'process_low'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
 
    } else {
        container "${params.yascp_container_docker}"
    }

    publishDir  path: "${outdir}/clustering_and_integration/plots",
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        val(outdir)
        path(file__filtered_cells)

    output:
        path("*.png") optional true
        path("*.pdf") optional true

    script:

        """
        echo "publish_directory: ${outdir}"
        plot_filtered_cells.py \
            --tsv_file ${file__filtered_cells} \
            --output_file adata-cell_filtered_per_experiment

        """
}



process PLOT_PCS {
    // Takes annData object with PCs and returns plots
    // ------------------------------------------------------------------------
    //cache false        // cache results from run
    scratch false      // use tmp directory

    label 'process_medium_memory'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
       container "${params.yascp_container_docker}"
    }

    publishDir  path: "${outdir}/clustering_and_integration/plots",
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        val(outdir)
        path(file__anndata)
        each n_pcs
        val(colors_quantitative)
        val(colors_categorical)

    output:
        val(outdir, emit: outdir)
        path("*.png"), emit: out_png
        path("*.pdf") optional true

    script:
        
        outfile = "pca"
        cmd__colors_quant = ""
        if (colors_quantitative != "") {
            cmd__colors_quant = "--colors_quantitative '${colors_quantitative}'"
        }
        cmd__colors_cat = ""
        if (colors_categorical != "") {
            cmd__colors_cat = "--colors_categorical '${colors_categorical}'"
        }

        """
        pca_plot.py \
            --h5_anndata ${file__anndata} \
            --num_pcs ${n_pcs} \
            ${cmd__colors_quant} \
            ${cmd__colors_cat} \
            --output_file ${outfile}
        """
}


process PLOT_PREDICTED_SEX {
    // Takes annData object, plots the predicted sex fron gene expression
    // ------------------------------------------------------------------------
    //cache false        // cache results from run
    label 'process_low'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }

    publishDir  path: "${outdir}/clustering_and_integration/plots",
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        val(outdir)
        path(file__anndata)

    output:
        path("*.png") optional true
        path("*.pdf") optional true

    script:

        outfile = "scatterplot-sex_sample_swap_check"
        """
        echo "publish_directory: ${outdir}"
        plot_predicted_sex.py \
            --h5_anndata ${file__anndata} \
            --output_file ${outfile}
        """
}

process PLOT_QC {
    // Takes annData object, generates basic qc plots
    // ------------------------------------------------------------------------
    //cache false        // cache results from run
    label 'process_low'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"

    }

    publishDir  path: "${outdir}/clustering_and_integration/plots",
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        val(outdir)
        path(file__anndata)
        each facet_columns

    output:
        path("*.png")
        path("*.pdf") optional true
        path("*.tsv") optional true

    script:

        outfile = "outfile"
        // Figure out if we are facetting the plot and update accordingly.
        cmd__facet_columns = ""
        if (facet_columns != "") {
            cmd__facet_columns = "--facet_columns ${facet_columns}"
        }
        """
        
            echo "publish_directory: ${outdir}"
            plot_qc_umi_nfeature_mt.py \
                --h5_anndata ${file__anndata} \
                --output_file ${outfile} \
                ${cmd__facet_columns}
            plot_qc_umi_mt_density.py \
                --h5_anndata ${file__anndata} \
                --output_file ${outfile} \
                ${cmd__facet_columns}
            plot_qc_nfeature_mt_density.py \
                --h5_anndata ${file__anndata} \
                --output_file ${outfile} \
                ${cmd__facet_columns}
            calculate_mads.py \
                --h5_anndata ${file__anndata} \
                --qc_key ${params.mads_categories} \
                --output_file mads
        """
}

process PLOT_DISTRIBUTIONS {
    // Takes annData object, generates basic qc plots
    // ------------------------------------------------------------------------
    //cache false        // cache results from run
    tag "${samplename}"
    label 'process_low'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }

    publishDir  path: "${outdir}/clustering_and_integration/plots",
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        val(outdir)
        path(file__anndata)
        each facet_columns
        each variable_columns_distribution_plots

    output:

        path("*.png")
        path("*.pdf") optional true
        path("*.tsv") optional true

    script:
        outfile = "outfile"
        // Figure out if we are facetting the plot and update accordingly.
        cmd__facet_columns = ""
        if (facet_columns != "") {
            cmd__facet_columns = "--facet_columns ${facet_columns}"
        }
        // Run distribution across cells if a value is specified
        cmd__anndataobs = ""
        cmd__anndataobs_ecdf = ""
        if (variable_columns_distribution_plots != "") {
            cmd__anndataobs = "plot_anndataobs_distribution_across_cells.py"
            cmd__anndataobs = "${cmd__anndataobs} --h5_anndata ${file__anndata}"
            cmd__anndataobs = "${cmd__anndataobs} --output_file ${outfile}"
            cmd__anndataobs = "${cmd__anndataobs} --variable_columns ${variable_columns_distribution_plots}"
            cmd__anndataobs = "${cmd__anndataobs} ${cmd__facet_columns}"
            cmd__anndataobs_ecdf = "${cmd__anndataobs} --ecdf"
        }
        """
        echo "publish_directory: ${outdir}"
        ${cmd__anndataobs}
        ${cmd__anndataobs_ecdf}
        """
}
