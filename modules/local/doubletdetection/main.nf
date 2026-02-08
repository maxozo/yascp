process DOUBLET_DETECTION {

    tag "${experiment_id}"
    label 'process_low'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }
    
    publishDir  path: "${params.outdir}/doublet_detection/multiplet.method=doubletdetection",
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        tuple(
            val(experiment_id),
            path(gex_h5ad)
        )

    output:
        path("plots/*.pdf") optional true
        path("plots/*.png") optional true
        tuple val(experiment_id), path("${experiment_id}__DoubletDetection_results.txt"), emit: result
        path "versions.yml", emit: versions

    script:
        
        outdir = "${params.outdir}/doublet_detection/multiplet"
        outdir = "${outdir}.method=doubletdetection"
        outfile = "${experiment_id}"

        """
            DoubletDetection.py --tenxdata_dir ${gex_h5ad} --n_iterations 100
            ln -s DoubletDetection_results.txt ${experiment_id}__DoubletDetection_results.txt

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                python: \$(python --version | sed 's/Python //g')
                python library argparse: \$(python -c "import argparse; print(argparse.__version__)")
                python library doubletdetection: \$(python -c "import doubletdetection; print(doubletdetection.__version__)")
                python library matplotlib: \$(python -c "import matplotlib; print(matplotlib.__version__)")
                python library numpy: \$(python -c "import numpy; print(numpy.__version__)")
                python library pandas: \$(python -c "import pandas; print(pandas.__version__)")
                python library scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
                python library tarfile: \$(python -c "import tarfile; print(tarfile.__version__)")
            END_VERSIONS
        """
}