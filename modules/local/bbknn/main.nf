// Default parameters for bbknn
params.bbknn = [
    run_process: false,
    batch_variable: [value: ["experiment_id"]]
]

def random_hex(n) {
    Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}


process BBKNN{
    
    // Calulates bbknn neighbors and saves UMAPS of these
    // ------------------------------------------------------------------------
    //cache false        // cache results from run
    scratch false      // use tmp directory

    // label 'process_medium'
    memory { 
            def sizeInGB = file__anndata.size() / 1e9 * 3 * task.attempt
            def minimumGB = 5
            return ((sizeInGB < minimumGB ? minimumGB : sizeInGB).toString() + 'GB')
        }

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
        each batch_var

    output:
        val(outdir, emit: outdir)
        path("${outfile}-bbknn.h5ad", emit: anndata)
        path(file__metadata, emit: metadata)
        path(file__pcs, emit: pcs)
        path("reduced_dims.tsv.gz", emit: reduced_dims)
        // NOTE: passing the param details as an unpublished file is messy,
        // but I could not get collect of ${param_details} and file to work.
        path(
            "reduced_dims-${param_details}.tsv.gz",
            emit: reduced_dims_params
        )
        path "versions.yml", emit: versions

    script:
        param_details = "${params__pcs}-bbknn"
        param_details = "${param_details}.batch=${batch_var}"
        param_details = "${param_details}.n_pcs=${n_pcs}"
        outdir = "${outdir_prev}/reduced_dims-${param_details}"
        // from the file__anndata job.
        outfile = "outfile_adata"
        """
        
        bbknn.py \
            --h5_anndata ${file__anndata} \
            --batch_key ${batch_var} \
            --n_pcs ${n_pcs} \
            --output_file ${outfile}-bbknn
        cp ${outfile}-bbknn-reduced_dims.tsv.gz \
            reduced_dims-${param_details}.tsv.gz
        mv ${outfile}-bbknn-reduced_dims.tsv.gz \
            reduced_dims.tsv.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/Python //g')
            python library argparse: \$(python -c "import argparse; print(argparse.__version__)")
            python library bbknn: \$(python -c "import bbknn; print(bbknn.__version__)")
            python library distutils: \$(python -c "import distutils; print(distutils.__version__)")
            python library numpy: \$(python -c "import numpy; print(numpy.__version__)")
            python library pandas: \$(python -c "import pandas; print(pandas.__version__)")
            python library scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
        END_VERSIONS
        """

}
