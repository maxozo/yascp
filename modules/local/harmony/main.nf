// Default parameters for harmony.
params.harmony = [
    run_process: false,
    variables_and_thetas: [value: [
        [variable: "experiment_id", theta: "1.0"],
        [variable: "experiment_id,phase", theta: "1.0,0.2"]
    ]]
]

def random_hex(n) {
    Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}

process HARMONY{
    // Takes PCs (rows = cell barcodes) and metadata (rows = cell barcodes),
    // runs Harmony
    // ------------------------------------------------------------------------
    //cache false        // cache results from run
    scratch false      // use tmp directory
    // storeDir '/tmp'

    memory { 
            def sizeInGB = file__anndata.size() / 1e9 * 3 * task.attempt
            def minimumGB = 5
            return ((sizeInGB < minimumGB ? minimumGB : sizeInGB).toString() + 'GB')
        }

    // label 'process_medium'
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
                        filename.replaceAll("${runid}-", "")
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
        each variables_and_thetas

    output:
        val(outdir, emit: outdir)
        path(file__anndata, emit: anndata)
        path(file__metadata, emit: metadata)
        path(file__pcs, emit: pcs)
        path("${runid}-reduced_dims.tsv.gz", emit: reduced_dims)
        // NOTE: passing the param details as an unpublished file is messy,
        // but I could not get collect of ${param_details} and file to work.
        path(
            "reduced_dims-${param_details}.tsv.gz",
            emit: reduced_dims_params
        )
        path "versions.yml", emit: versions


    script:
        runid = "${params__pcs}_${n_pcs}_${variables_and_thetas.variable}"
        param_details = "${params__pcs}-harmony"
        param_details = "${param_details}.n_pcs=${n_pcs}"
        param_details = "${param_details}.variables=${variables_and_thetas.variable}"
        outdir = "${outdir_prev}/reduced_dims-${param_details}"
        outdir = "${outdir}.thetas=${variables_and_thetas.theta}"
        theta_str = "${variables_and_thetas.theta}".replaceAll("\\.", "pt")
        param_details = "${param_details}.thetas=${theta_str}"

        """
            harmony_process_pcs.py \
                --pca_file ${file__pcs} \
                --metadata_file ${file__metadata} \
                --metadata_columns ${variables_and_thetas.variable} \
                --theta ${variables_and_thetas.theta} \
                --n_pcs ${n_pcs} \
                --out_file ${runid}-reduced_dims
            cp ${runid}-reduced_dims.tsv.gz \
            reduced_dims-${param_details}.tsv.gz

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                python: \$(python --version | sed 's/Python //g')
                python library argparse: \$(python -c "import argparse; print(argparse.__version__)")
                python library csv: \$(python -c "import csv; print(csv.__version__)")
                python library distutils: \$(python -c "import distutils; print(distutils.__version__)")
                python library harmonypy: \$(python -c "import harmonypy; print(harmonypy.__version__)")
                python library numpy: \$(python -c "import numpy; print(numpy.__version__)")
                python library pandas: \$(python -c "import pandas; print(pandas.__version__)")
            END_VERSIONS
        """
        // NOTE: below code for harmony in R
        // harmony_process_pcs.R \
        //     --pca_file ${file__pcs} \
        //     --metadata_file ${file__metadata} \
        //     --metadata_columns ${variables_and_thetas.variable} \
        //     --theta ${variables_and_thetas.theta} \
        //     --n_pcs ${n_pcs} \
        //     --out_file reduced_dims
}
