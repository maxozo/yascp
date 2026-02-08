def random_hex(n) {
    Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}


if (binding.hasVariable("echo_mode") == false) {
    echo_mode = true
}


process CLUSTER {
    // Clusters results.
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
                    if (filename.endsWith("normalized_pca.h5ad")) {
                        null
                    } else if(filename.endsWith("metadata.tsv.gz")) {
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
        val(outdir_prev)
        path(file__anndata)
        path(file__metadata)
        path(file__pcs)
        path(file__reduced_dims)
        each number_neighbors
        each method
        each resolution
        // tuple(val(outdir_prev), path(file__reduced_dims))

    output:
        val(outdir, emit: outdir)
        path("${outfile}-clustered.h5ad", emit: anndata)
        path(file__metadata, emit: metadata)
        path(file__pcs, emit: pcs)
        path(file__reduced_dims, emit: reduced_dims)
        path("${outfile}-clustered.tsv.gz", emit: clusters)
        val(outdir_prev, emit: outdir__reduced_dims)
        path("*.pdf") optional true
        path("*.png") optional true
        
    script:
        
        resolution_str = "${resolution}" //.replaceAll("\\.", "pt")
        outdir = "${outdir_prev}/cluster"
        outdir = "${outdir}.number_neighbors=${number_neighbors}"
        outdir = "${outdir}.method=${method}"
        outdir = "${outdir}.resolution=${resolution_str}"
        // from the file__anndata job.
        outfile = "clustering_${resolution}"

        """
        echo "publish_directory: ${outdir}"
        scanpy_cluster.py \
            --h5_anndata ${file__anndata} \
            --tsv_pcs ${file__reduced_dims} \
            --number_neighbors ${number_neighbors} \
            --cluster_method ${method} \
            --resolution ${resolution} \
            --number_cpu ${task.cpus} \
            --output_file ${outfile}-clustered
        """
}

process PLOT_PHENOTYPE_ACROSS_CLUSTERS {
    // Takes annData object, plots distribution of obs value across clusters
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
                saveAs: {filename -> filename.replaceAll("-", "")},
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(file__anndata)
        each variables

    output:
        val(outdir, emit: outdir)
        path("plots/*.png")
        path("plots/*.pdf") optional true

    script:
        
        outdir = "${outdir_prev}"
        // from the file__anndata job.
        outfile = "${file__anndata}".minus(".h5ad")
            .split("-").drop(1).join("-")
        // Append run_id to output file.
        outfile = "${variables}-adata-cluster_boxplot"
        """

        rm -fr plots
        plot_anndataobs_across_clusters.py \
            --h5_anndata ${file__anndata} \
            --pheno_columns ${variables} \
            --output_file ${outfile}
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}

process SERIALIZE_KNOWN_MARKERS {
    // Serializes known markers for analysis
    // ------------------------------------------------------------------------
    scratch false      // use tmp directory
    label 'process_low'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"

    } else {
        container "${params.yascp_container_docker}"
    }

    input:
        tuple(
            val(file_id),
            file(marker_file)
        )

    output:
        path("${outfile}", emit: marker_file)

    script:
        
        outfile = "${file_id}"
        cmd__run = "touch NO_FILE.tsv"
        if (file_id != "") {
            outfile = "${file_id}.tsv"
            cmd__run = "ln --physical ${marker_file} ${outfile}"
        } else {
            outfile = "NO_FILE"
            cmd__run = "touch ${outfile}"
        }
        """
        ${cmd__run}
        """
}

process PLOT_KNOWN_MARKERS {
    // Plots markers from previous studies as dotplots
    // ------------------------------------------------------------------------
    //cache false        // cache results from run
    scratch false        // use tmp directory
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
        each file(marker_file)

    output:
        val(outdir, emit: outdir)
        path("plots_known_markers/*.pdf") optional true
        path("plots_known_markers/*.png") optional true

    script:
        
        outdir = "${outdir_prev}"
        // from the file__anndata job.
        // outfile = "${file__anndata}".minus(".h5ad")
        //     .split("-").drop(1).join("-")
        outfile = "${marker_file.name}".minus(".tsv")
        // Only run this script if there is a value
        cmd__run = ""
        if (outfile != "NO_FILE") {
            cmd__run = "plot_known_markers.py"
            cmd__run = "${cmd__run} --h5_anndata ${file__anndata}"
            cmd__run = "${cmd__run} --markers_database ${marker_file}"
            cmd__run = "${cmd__run} --output_file ${outfile}"
        }
        """
        rm -fr plots_known_markers
        ${cmd__run}
        mkdir plots_known_markers
        mv *pdf plots_known_markers/ 2>/dev/null || true
        mv *png plots_known_markers/ 2>/dev/null || true
        """
}

process CLUSTER_VALIDATE_RESOLUTION_KERAS {
    // Validate the resolution for clusters.
    // ------------------------------------------------------------------------
    //cache false        // cache results from run
    //maxForks 2         // hard to control memory usage. limit to 3 concurrent
    label 'process_low'        // use GPU
    scratch false      // use tmp directory
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
        
    } else {
        container "${params.yascp_container_docker}"
    }

    publishDir  path: "${outdir}",
                saveAs: { filename ->
                    if (filename.endsWith("clustered.h5ad") ||
                        filename.endsWith("metadata.tsv.gz") ||
                        filename.endsWith("pcs.tsv.gz") ||
                        filename.endsWith("reduced_dims.tsv.gz") ||
                        filename.endsWith("clustered.tsv.gz")) {
                        return null
                    }

                    def parts = filename.split("___")
                    return (parts.length > 1) ? parts[1] : filename
                },
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(file__anndata)
        path(file__metadata)
        path(file__pcs)
        path(file__reduced_dims)
        path(file__clusters)
        each sparsity
        each train_size_cells
        val(outdir__reduced_dims)

    output:
        val(outdir, emit: outdir)
        path(file__anndata, emit: anndata)
        path(file__metadata, emit: metadata)
        path(file__pcs, emit: pcs)
        path(file__reduced_dims, emit: reduced_dims)
        path(file__clusters, emit: clusters)
        path("*${outfile}.h5", emit: model)
        path("*${outfile}.json", emit: model_yaml)
        path("*weights.h5", emit: model_weights)
        path("*${outfile}-model_report.tsv.gz", emit: model_report)
        path(
            "*${outfile}-test_result.tsv.gz",
            emit: model_test_result
        )
        path(
            "*${outfile}-weights.tsv.gz",
            emit: model_weights_tsv
        )
        tuple(
            val("${outdir__reduced_dims}"),
            file("*${outfile}-model_report.tsv.gz",),
            file("*${outfile}-test_result.tsv.gz",),
            emit: plot_input
        )
        path("plots/*.png") optional true
        path("plots/*.pdf") optional true
        path "versions.yml", emit: versions

    script:
        outdir = "${outdir_prev}/validate_resolution"
        // outdir = "${outdir}.method=${method}"
        // from the file__anndata job.
        outfile = "${file__anndata}".minus(".h5ad")
            .split("-").drop(1).join("-")
        // Add sparsity information.
        sparsity_txt = "${sparsity}".replaceAll("\\.", "pt")
        outfile = "${outfile}-sparsity_l1=${sparsity_txt}"
        // Add info on number of cells for training
        // cmd__train_cells = ""
        // if (train_size_cells > 0) {
        //     train_size_cells_txt = "${train_size_cells}"
        //     cmd__train_cells = "--train_size_cells ${train_size_cells}"
        // }
        outfile = "${outfile}-train_size_cells=${train_size_cells}"

        """
        vcf_name=\$(python ${projectDir}/bin/random_id.py)
        rm -fr plots
        scanpy_cluster_validate_resolution-keras.py \
            --h5_anndata ${file__anndata} \
            --sparsity_l1 ${sparsity} \
            --number_epoch 25 \
            --batch_size 32 \
            --train_size_cells ${train_size_cells} \
            --output_file \${vcf_name}___${outfile}
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/Python //g')
            python library argparse: \$(python -c "import argparse; print(argparse.__version__)")
            python library csv: \$(python -c "import csv; print(csv.__version__)")
            python library distutils: \$(python -c "import distutils; print(distutils.__version__)")
            python library keras: \$(python -c "import keras; print(keras.__version__)")
            python library matplotlib: \$(python -c "import matplotlib; print(matplotlib.__version__)")
            python library numpy: \$(python -c "import numpy; print(numpy.__version__)")
            python library pandas: \$(python -c "import pandas; print(pandas.__version__)")
            python library plotnine: \$(python -c "import plotnine; print(plotnine.__version__)")
            python library scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
            python library scikeras: \$(python -c "import scikeras; print(scikeras.__version__)")
            python library scipy: \$(python -c "import scipy; print(scipy.__version__)")
            python library sklearn: \$(python -c "import sklearn; print(sklearn.__version__)")
            python library tensorflow: \$(python -c "import tensorflow; print(tensorflow.__version__)")
        END_VERSIONS
        """
}

process PLOT_RESOLUTION_VALIDATE {
    // Plot the AUC from validation models across resolutions
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
                    if (filename.endsWith("clustered.h5ad")) {
                        null
                    } else if(filename.endsWith("metadata.tsv.gz")) {
                        null
                    } else if(filename.endsWith("pcs.tsv.gz")) {
                        null
                    } else if(filename.endsWith("reduced_dims.tsv.gz")) {
                        null
                    } else if(filename.endsWith("clustered.tsv.gz")) {
                        null
                    } else {
                        filename.replaceAll("-", "")
                    }
                },
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        tuple(
            val(outdir_prev),
            path(files__model_report),
            path(files__y_prob_df)
        )

    output:
        val(outdir, emit: outdir)
        path(
            "${outfile}-merged_model_report.tsv.gz",
            emit: merged_model_report
        )
        path(
            "${outfile}-merged_test_result.tsv.gz",
            emit: merged_test_result
        )
        path("plots/*.png") optional true
        path("plots/*.pdf") optional true

    script:
        outdir = "${outdir_prev}"
        files__model_report = files__model_report.join('::')
        files__y_prob_df = files__y_prob_df.join('::')
        outfile = "resolution_tuning"
        """
        echo "publish_directory: ${outdir}"
        rm -fr plots
        plot_resolution_boxplot.py \
            --model_reports ${files__model_report} \
            --h_line 0.75 \
            --output_file ${outfile}-hline0pt8
        plot_resolution_boxplot.py \
            --model_reports ${files__model_report} \
            --output_file ${outfile}
        plot_resolution_curve.py \
             --y_prob_dfs ${files__y_prob_df} \
             --output_file ${outfile}
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}

process CLUSTER_MARKERS {
    // Find markers for clusters.
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
                    if (filename.endsWith("clustered.h5ad")) {
                        null
                    } else if(filename.endsWith("metadata.tsv.gz")) {
                        null
                    } else if(filename.endsWith("pcs.tsv.gz")) {
                        null
                    } else if(filename.endsWith("reduced_dims.tsv.gz")) {
                        null
                    } else if(filename.endsWith("clustered.tsv.gz")) {
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
        path(file__reduced_dims)
        path(file__clusters)
        each method

    output:
        val(outdir, emit: outdir)
        path(file__anndata, emit: anndata)
        path(file__metadata, emit: metadata)
        path(file__pcs, emit: pcs)
        path(file__reduced_dims, emit: reduced_dims)
        path(file__clusters, emit: clusters)
        path(
            "${outfile}-cluster_markers.tsv.gz",
            emit: cluster_markers
        )
        path(
            "${outfile}-cluster_markers-filter__*.tsv.gz"
        ) optional true
        path("plots/*.pdf") optional true
        path("plots/*.png") optional true

    script:
        outdir = "${outdir_prev}/cluster_markers"
        outdir = "${outdir}.method=${method}"
        // from the file__anndata job.
        outfile = "${file__anndata}".minus(".h5ad")
            .split("-").drop(1).join("-")
        // Only run this script if there is a value
        cmd__enrich = ""
        // if (method == "wilcoxon") {
        //     f = "${outfile}-cluster_markers"
        //     cmd__enrich = "0056-marker_enrichmet.R"
        //     cmd__enrich = "${cmd__enrich} --markers_table ${f}.tsv.gz;"
        //     cmd__enrich = "${cmd__enrich} 0056-marker_enrichmet.R"
        //     cmd__enrich = "${cmd__enrich} --markers_table ${f}-filter__fdr0pt05.tsv.gz;"
        //     cmd__enrich = "${cmd__enrich} 0056-marker_enrichmet.R"
        //     cmd__enrich = "${cmd__enrich} --markers_table ${f}-filter__fdr0pt05__user_variable_genes_exclude.tsv.gz;"
        // }

        """
        rm -fr plots
        scanpy_cluster_markers.py \
            --h5_anndata ${file__anndata} \
            --rank_genes_method ${method} \
            --number_cpu ${task.cpus} \
            --output_file ${outfile}
        ${cmd__enrich}
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}

process CELLEX_CLUSTER_MARKERS {
    // Find markers for clusters using CELLEX.
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
                saveAs: {filename -> filename.replaceAll("-", "")},
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(file__anndata)

    output:
        path(
            "${outfile}-essd.tsv.gz"
        )
        path(
            "${outfile}-esmu*.tsv.gz"
        )
    when:
        params.cellex_cluster_markers

        
    script:
    // The script generates two files, adata-normalized_pca-bbknn-umap-clustered.esmu.csv.gz
    // and adata-normalized_pca-bbknn-umap-clustered.essd.csv.gz.
    // These are two matrices in the format genes x clusters. The first file, esmu.csv.gz,
    // contains the expression specificity scores for each gene within each cluster (scale from 0 to 1),
    // while the second file essd.csv.gz comprises the respective standard deviation for each gene and cluster.
        outdir = "${outdir_prev}/cluster_markers.method=cellex"
        // from the file__anndata job.
        outfile = "${file__anndata}".minus(".h5ad")
            .split("-").drop(1).join("-")
        // Only run this script if there is a value

        """

        cellex_cluster_markers.py \
            --h5_anndata ${file__anndata} \
            --output_file ${outfile} \
            --verbose False
        """
}

process PREP_CELLXGENE {
    // Preps adata file for cellxgene
    // ------------------------------------------------------------------------
    //cache false           // cache results from run
    scratch false           // use tmp directory
    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"

    } else {
        container "${params.yascp_container_docker}"
    }

    publishDir  path: "${outdir}",
                saveAs: {filename ->
                    if (filename.endsWith("clustered.h5ad")) {
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

    output:
        val(outdir, emit: outdir)
        path(
            "${outfile}-cellxgene.h5ad",
            emit: cluster_markers
        )
        path "versions.yml", emit: versions

    script:
        outdir = "${outdir_prev}"
        // from the file__anndata job.
        outfile = "${file__anndata}".minus(".h5ad").split("-").drop(1).join("-")

        """

        cellxgene.py \
            --h5_anndata ${file__anndata} \
            --drop_extra_info \
            --output_file ${outfile}-cellxgene

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/Python //g')
            python library argparse: \$(python -c "import argparse; print(argparse.__version__)")
            python library numpy: \$(python -c "import numpy; print(numpy.__version__)")
            python library scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
            python library scipy: \$(python -c "import scipy; print(scipy.__version__)")
            python library warnings: \$(python -c "import warnings; print(warnings.__version__)")
        END_VERSIONS
        """
}

process CONVERT_SEURAT {
    // Converts anndata h5 file to a Seurat data object.
    // TODO: automatically add reduced_dims to Seurat data object.
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
                saveAs: {filename -> filename.replaceAll("-", "")},
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(file__anndata)

    output:
        val(outdir, emit: outdir)
        path("${outfile}.rds.gz", emit: seurat_data)
        path("${outdir_relative}/*", emit: matrix_data)
        // tuple(path("${outdir_relative}/*"), emit: matrix_data)

    script:
        // from the file__anndata job.
        outfile = "${file__anndata}".minus(".h5ad").split("-").drop(1).join("-")
        outdir_relative = "matrices-${outfile}"
        outdir = "${outdir_prev}"

        """
        echo "publish_directory: ${outdir}"
        outdir_relative_full_path=\$(pwd)"/${outdir_relative}"
        convert-anndata_10x.py \
            --h5_anndata ${file__anndata} \
            --output_dir ${outdir_relative}
        ln -s \${outdir_relative_full_path}/matrix-X.mtx.gz ${outdir_relative}/matrix.mtx.gz
        convert-10x_seurat.R \
            --in_dir \${outdir_relative_full_path} \
            --metadata_file \${outdir_relative_full_path}/metadata-barcodes.tsv.gz \
            --count_matrix_file \${outdir_relative_full_path}/matrix-counts.mtx.gz \
            --out_file ${outfile}
        unlink ${outdir_relative}/matrix.mtx.gz
        """
}
