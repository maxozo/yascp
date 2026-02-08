

def random_hex(n) {
    Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}

if (binding.hasVariable("echo_mode") == false) {
    echo_mode = true
}

process SCRUBLET {
    // Runs scrublet for each sample.

    tag "${experiment_id}"

    label 'process_low'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }
    
    publishDir  path: "${params.outdir}/doublet_detection/scrublet",
                saveAs: {filename ->
                    if (filename.endsWith("multiplet_calls_published.txt")) {
                        null
                    }
                    else if (filename.endsWith(".gz")){
                        null
                    }
                    else {
                        filename.replaceAll("-", "")
                    }
                },
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        tuple(
            val(experiment_id),
            path(file_10x)
        )
        val(expected_multiplet_rate)
        val(n_simulated_multiplet)
        val(multiplet_threshold_method)
        val(scale_log10)

    output:
        // val(outdir, emit: outdir)
        tuple val(experiment_id), path("${outfile}-scrublet.tsv.gz"), emit: scrublet_paths
        val(experiment_id, emit: experiment_id)
        path("${outfile}-scrublet.tsv.gz", emit: multiplet_calls)
        path(
            "${outfile}-multiplet_calls_published.txt",
            emit: multiplet_calls_published
        )
        path("plots/*.pdf") optional true
        path("plots/*.png") optional true
        tuple val(experiment_id), path("${experiment_id}-scrublet.tsv"), emit: result
        path "versions.yml", emit: versions

    script:
        
        outdir = "${params.outdir}/doublet_detection/"
        outdir = "${outdir}scrublet"
        outfile = "${experiment_id}"
        // Check to see if we should use use log10 of the doublet simulations
        // to derive the threshold
        cmd__scale_log10 = ""
        if (scale_log10 == "True") {
            cmd__scale_log10 = "--scale_log10"
        }

        """

        rm -fr plots
        run_scrublet.py \
            --tenxdata_dir ${file_10x} \
            --expected_multiplet_rate ${expected_multiplet_rate} \
            --n_simulated_multiplet ${n_simulated_multiplet} \
            --multiplet_threshold_method ${multiplet_threshold_method} \
            ${cmd__scale_log10} \
            --output_file ${outfile}
            
        echo -e "${experiment_id}\t${outdir}/${outfile}-scrublet.tsv.gz" > \
            ${outfile}-multiplet_calls_published.txt
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true

        zcat ${experiment_id}-scrublet.tsv.gz | awk -F' ' '{print \$1"\\t"\$3"\\t"\$2}' > ${experiment_id}-scrublet.tsv

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/Python //g')
            python library argparse: \$(python -c "import argparse; print(argparse.__version__)")
            python library csv: \$(python -c "import csv; print(csv.__version__)")
            python library distutils: \$(python -c "import distutils; print(distutils.__version__)")
            python library matplotlib: \$(python -c "import matplotlib; print(matplotlib.__version__)")
            python library numpy: \$(python -c "import numpy; print(numpy.__version__)")
            python library pandas: \$(python -c "import pandas; print(pandas.__version__)")
            python library plotnine: \$(python -c "import plotnine; print(plotnine.__version__)")
            python library scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
            python library scrublet: \$(pip show scrublet | grep '^Version:' | cut -d':' -f2)
            python library seaborn: \$(python -c "import seaborn; print(seaborn.__version__)")
            python library skimage: \$(python -c "import skimage; print(skimage.__version__)")
        END_VERSIONS
        """
}
