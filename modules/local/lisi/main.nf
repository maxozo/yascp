
def random_hex(n) {
    Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}


process LISI{
    // Takes a list of reduced_dims and calculates lisi
    // ------------------------------------------------------------------------
    //cache false        // cache results from run
    scratch false      // use tmp directory
    // label 'process_high'
    label 'process_medium'
    
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
                    } else {
                        filename.replaceAll("-", "")
                    }
                },
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(file__metadata)
        val(variables)
        file(file__reduced_dims)
        //tuple(val(label__reduced_dims), file(file__reduced_dims))

    output:
        val(outdir, emit: outdir)
        path(file__metadata, emit: metadata)
        path("${outfile}-lisi.tsv.gz", emit: clusters)
        path("plots/*.pdf") optional true
        path("plots/*.png") optional true
        path "versions.yml", emit: versions

    script:
        outdir = "${outdir_prev}"

        // from the file__metadata job.
        outfile = "${file__metadata}".minus(".tsv.gz")
            .split("-").drop(1).join("-")
        file__reduced_dims = file__reduced_dims.join("::")
        label__reduced_dims = file__reduced_dims
            .replaceAll("reduced_dims-", "")
            .replaceAll(".tsv.gz", "")

        """
        rm -fr plots
        lisi.py \
            --reduced_dims_tsv ${file__reduced_dims} \
            --reduced_dims_tsv_labels ${label__reduced_dims} \
            --metadata_tsv ${file__metadata} \
            --metadata_columns ${variables} \
            --perplexity 30 \
            --output_file ${outfile}-lisi
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/Python //g')
            python library argparse: \$(python -c "import argparse; print(argparse.__version__)")
            python library csv: \$(python -c "import csv; print(csv.__version__)")
            python library distutils: \$(python -c "import distutils; print(distutils.__version__)")
            python library harmonypy: \$(python -c "import harmonypy; print(harmonypy.__version__)")
            python library pandas: \$(python -c "import pandas; print(pandas.__version__)")
            python library plotnine: \$(python -c "import plotnine; print(plotnine.__version__)")
        END_VERSIONS
        """

}
