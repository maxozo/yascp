
def random_hex(n) {
    Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}

if (binding.hasVariable("echo_mode") == false) {
    echo_mode = true
}

process ESTIMATE_PCA_ELBOW {
    // Takes annData object, estiamtes the elbow in PC var explained.
    // ------------------------------------------------------------------------
    //cache false        // cache results from run
    scratch false      // use tmp directory
    echo echo_mode          // echo output from script

    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"

    } else {
        container "${params.yascp_container_docker}"
    }

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("-", "")},
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(file__anndata)
        val(add_n_to_estimate)

    output:
        val(outdir, emit: outdir)
        path("${outfile}.tsv", emit: pca_elbow_estimate)
        env(AUTO_ELBOW, emit: auto_elbow)
        path("plots/*.png")
        path("plots/*.pdf") optional true

    script:
        
        outdir = "${outdir_prev}"
        outfile = "${file__anndata}".minus(".h5ad")
        outfile = "${outfile}-knee"
        """
            rm -fr plots
            estimate_pca_elbow.py \
                --h5_anndata ${file__anndata} \
                --add_n_pcs_to_elbow ${add_n_to_estimate} \
                --output_file ${outfile}
            mkdir plots
            mv *pdf plots/ 2>/dev/null || true
            mv *png plots/ 2>/dev/null || true
            AUTO_ELBOW=\$(cat ${outfile}-auto_elbow_estimate.tsv)
        """
}
