#!/usr/bin/env nextflow

def random_hex(n) {
    Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}

if (binding.hasVariable("echo_mode") == false) {
    echo_mode = true
}

process UMAP_CALCULATE {
    // UMAP from reduced_dims.
    // ------------------------------------------------------------------------
    //cache false        // cache results from run
    scratch false      // use tmp directory
    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }

    input:
        val(outdir)
        path(file__anndata)
        path(file__metadata) // not really needed
        path(file__pcs) // not really needed
        path(file__reduced_dims)
        val(use_pcs_as_reduced_dims) // For BBKNN
        each n_neighbors
        each umap_init
        each umap_min_dist
        each umap_spread
        val(method)

    output:
        val(outdir, emit: outdir)
        tuple(
            val("${in_file_id}"),
            val("${outdir}"),
            file(file__anndata),
            file(file__metadata),
            file(file__pcs),
            file(file__reduced_dims),
            file("*${outfile_pattern}.h5ad"),
            emit: outdir_anndata
        )
        path("*${outfile_pattern}.h5ad",emit:adata_out)

    script:
        runid = random_hex(16)
        in_file_id = file__reduced_dims.toString().tokenize('-').get(0)
        // For output file, use anndata name. First need to drop the runid
        // from the file__anndata job.
        outfile = "${file__anndata}".minus(".h5ad")
        outfile_pattern = "${outfile}_${method}_umap_${n_neighbors}_${umap_spread}_${umap_init}_${umap_min_dist}"
        
        outfile = "${runid}-${outfile}_${method}_umap_${n_neighbors}_${umap_spread}_${umap_init}_${umap_min_dist}"
        // Check to see if we should use PCs in the reduced dims slot.
        // Important for BBKNN where reduced_dims == UMAPs not adjusted PCs.
        cmd__tsv_pcs = "--tsv_pcs ${file__reduced_dims}"
        if (use_pcs_as_reduced_dims == "True") {
            cmd__tsv_pcs = "--tsv_pcs ${file__pcs}"
        }

        """
            echo ${in_file_id}
            umap_calculate.py \
                --h5_anndata ${file__anndata} \
                ${cmd__tsv_pcs} \
                --n_neighbors ${n_neighbors} \
                --umap_init ${umap_init} \
                --umap_min_dist ${umap_min_dist} \
                --umap_spread ${umap_spread} \
                --number_cpu ${task.cpus} \
                --output_file ${outfile}
        """
        //--calculate_densities \
}

process GENERATE_FINAL_UMAPS{

    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }

  publishDir  path: "${outdir}",
              mode: "${params.copy_mode}",
              overwrite: "true"
  input:
    path(file__anndata)
    val(outdir_prev)
  output:
    path("umap-*")
  
  script:
    outdir = "${outdir_prev}/handover/UMAPs"
    """
        umap_plot_final.py \
            --h5_anndata ${file__anndata} \
            --number_cpu 1 \
            --colors_quantitative '${params.umap.colors_quantitative.value}' \
            --colors_categorical '${params.umap.colors_categorical.value}' \
            --drop_legend_n 40 \
            --output_file UMAP
    """

}

process UMAP_GATHER {
    // Merge UMAP from reduced_dims (reduce or gather).
    // ------------------------------------------------------------------------
    //cache false        // cache results from run
    scratch false      // use tmp directory
    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }
    publishDir  path: "${outdir}",
                saveAs: {filename ->
                    if(filename.endsWith("metadata.tsv.gz")) {
                        null
                    } else if(filename.endsWith("pcs.tsv.gz")) {
                        null
                    } else if(filename.endsWith("reduced_dims.tsv.gz")) {
                        null
                    } else {
                        filename.replaceAll("-", "")
                    }
                },
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        tuple(
            val(key),
            val(outdir),
            path(original__file__anndata),
            path(original__file__metadata),
            path(original__file__pcs),
            path(original__file__reduced_dims),
            path(files__anndata)
        )
        
    output:
        val(outdir, emit: outdir)
        path("${outfile}.h5ad", emit: anndata)
        path(original__file__metadata, emit: metadata)
        path(original__file__pcs, emit: pcs)
        path(original__file__reduced_dims, emit: reduced_dims)

    script:
        
        outfile = "umap_gather_out"
        files__anndata = files__anndata.join(',')

        """
            umap_gather.py \
                --h5_root ${original__file__anndata} \
                --output_file ${outfile} \
                --h5_anndata_list ${files__anndata}
        """
}


process UMAP_PLOT_SWARM {
    // Plot UMAPs.
    // ------------------------------------------------------------------------
    //cache false        // cache results from run
    scratch false      // use tmp directory
    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }
    publishDir  path: "${outdir}/plots",
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        val(outdir)
        path(file__anndata)
        val(colors_quantitative)
        val(colors_categorical)
        val(drop_legend_n)

    output:
        path("*.png")
        path("*.pdf") optional true
    script:

        outfile = "umap"
        cmd__colors_quant = ""
        if (colors_quantitative != "") {
            cmd__colors_quant = "--colors_quantitative '${colors_quantitative}'"
        }
        cmd__colors_cat = ""
        if (colors_categorical != "") {
            cmd__colors_cat = "--colors_categorical '${colors_categorical}'"
        }

        """
            umap_plot.py \
                --h5_anndata ${file__anndata} \
                --number_cpu ${task.cpus} \
                ${cmd__colors_quant} \
                ${cmd__colors_cat} \
                --drop_legend_n ${drop_legend_n} \
                --output_file ${outfile}

        """
}


process UMAP_CALCULATE_AND_PLOT {
    // UMAP from reduced_dims.
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }
    publishDir  path: "${outdir}/plots",
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        val(outdir)
        path(file__anndata)
        path(file__pcs)
        path(file__reduced_dims)
        val(use_pcs_as_reduced_dims) // For BBKNN
        val(colors_quantitative)
        val(colors_categorical)
        each n_neighbors
        each umap_init
        each umap_min_dist
        each umap_spread

    output:
        path("*.png"), emit: dummy_output
        path("*.pdf") optional true

    script:
        outfile = "umap"
        cmd__colors_quant = ""
        if (colors_quantitative != "") {
            cmd__colors_quant = "--colors_quantitative '${colors_quantitative}'"
        }
        cmd__colors_cat = ""
        if (colors_categorical != "") {
            cmd__colors_cat = "--colors_categorical '${colors_categorical}'"
        }
        drop_legend_n = "-1"
        if (cmd__colors_cat.contains("experiment_id")) {
            drop_legend_n = "40"
        }
        cmd__tsv_pcs = "--tsv_pcs ${file__reduced_dims}"
        if (use_pcs_as_reduced_dims == "True") {
            cmd__tsv_pcs = "--tsv_pcs ${file__pcs}"
        }

        """

        umap_calculate_and_plot.py \
            --h5_anndata ${file__anndata} \
            ${cmd__tsv_pcs} \
            --n_neighbors ${n_neighbors} \
            --umap_init ${umap_init} \
            --umap_min_dist ${umap_min_dist} \
            --umap_spread ${umap_spread} \
            --number_cpu ${task.cpus} \
            ${cmd__colors_quant} \
            ${cmd__colors_cat} \
            --drop_legend_n ${drop_legend_n} \
            --output_file ${outfile}

        """
}
