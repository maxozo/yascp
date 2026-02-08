def random_hex(n) {
    Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}

process SUBSET_PCS{
   
    // Takes PCs (rows = cell barcodes) and subsets down to a specified number.
    // ------------------------------------------------------------------------
    //cache false        // cache results from run
    scratch false      // use tmp directory
    label 'process_low'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }


    publishDir  path: "${outdir}",
                saveAs: {filename ->
                    if (filename.endsWith("normalized_pca.h5ad")) {
                        null
                    } else if(filename.endsWith("metadata.tsv.gz")) {
                        null
                    } else if(filename.endsWith("pcs.tsv.gz")) {
                        null
                    } else if(filename.endsWith("${param_details}.tsv.gz")) {
                        null
                    } else {
                        filename.replaceAll("-", "")
                    }
                },
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(file__anndata)
        path(file__metadata)
        path(file__pcs)
        val(params__pcs)
        each n_pcs

    output:
        val(outdir, emit: outdir)
        path(file__anndata, emit: anndata)
        path(file__metadata, emit: metadata)
        path(file__pcs, emit: pcs)
        path("reduced_dims.tsv.gz", emit: reduced_dims)
        // NOTE: passing the param details as an unpublished file is messy,
        // but I could not get collect of ${param_details} and file to work.
        path(
            "reduced_dims-${param_details}.tsv.gz",
            emit: reduced_dims_params
        )

    script:
        
        param_details = "${params__pcs}-pca"
        param_details = "${param_details}.n_pcs=${n_pcs}"
        outdir = "${outdir_prev}/reduced_dims-${param_details}"
        """
        echo "publish_directory: ${outdir}"
        subset_pca_file.py \
            --tsv_pcs ${file__pcs} \
            --number_pcs ${n_pcs} \
            --output_file reduced_dims
        cp reduced_dims.tsv.gz \
            reduced_dims-${param_details}.tsv.gz
        """

}
