params {

    lisi{
        run_process=true
    }
    replace_genotype_ids=false
    cluster_validate_resolution_keras = true
    filter_outliers = false
    extra_sample_metadata =""
    output_dir = outdir= "${launchDir}/recluster_resolutions"
    cellex_cluster_markers=true 
    cluster_markers = false
    normalise_andata = false
    skip_handover = true
    split_ad_per_bach=true //if not splitting the celltype assignment will be run on full tranche
    skip_preprocessing = true

    harmony{
        run_process= true
    }
    umap{
        run_process = true
        colors_quantitative{
            description = 'Comma separated string of quantitative variables that will be used to color points.'
            value = 'n_cells,total_counts,pct_counts_gene_group__mito_transcript,prob_doublet,pct_counts_gene_group__ribo_rna,Azimuth:predicted.celltype.l2.score,Azimuth:mapping.score,log10_ngenes_by_count'
        }
        colors_categorical{
            description = 'Comma separated string of categorical variables that will be used to color points.'
            value = 'cell_passes_qc,cell_passes_qc-per:Azimuth:L0_predicted.celltype.l2,experiment_id,Azimuth:predicted.celltype.l2,Celltypist:Immune_All_Low:predicted_labels,Celltypist:Immune_All_High:predicted_labels,donor_id'
        }
    }

    mads_categories ='pct_counts_gene_group__mito_transcript,pct_counts_gene_group__mito_protein,pct_counts_gene_group__ribo_protein,pct_counts_gene_group__ribo_rna,total_counts,n_genes_by_counts,log10_ngenes_by_count'

    cluster{
        description = """Parameters for clustering. All pairwise combinations of
        method and resolution will be performed."""
        number_neighbors{
            description = """Number of neighbors. If <= 0, uses number of unique
            experiment_id."""
            value = 15
        }
        methods{
            description = 'Clustering method. Valid options [leiden|louvain].'
            value = 'leiden'
        }
        resolutions{
            description = 'Clustering resolution.'
            value = [0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
        }

        variables_boxplot{
            decription = 'Generate boxplots of these variables for each cluster.'
            value ='n_cells,total_counts,pct_counts_gene_group__mito_transcript'
        }

        known_markers{
            run_process = false
            description = """Files with markers that will be used to generate
            dotplots. Each marker file should be the full path and have the
            following columns: cell_type, hgnc_symbol. The following columns
            are optional: p_value_adj. Use "" for a single entry in the
            file_id and file value to indicate no plots."""
            value = [
                [ file_id: 'SmillieCS_31348891', file: '/lustre/scratch119/humgen/projects/sc-eqtl-ibd/data/marker_gene_db-raw_data/database/celltypes/colon/SmillieCS-31348891/database.tsv' ],
                [ file_id: 'ParikhK_30814735', file: '/lustre/scratch119/humgen/projects/sc-eqtl-ibd/data/marker_gene_db-raw_data/database/celltypes/colon/ParikhK-30814735/database.tsv' ],
                [ file_id: 'JamesKR_32066951', file: '/lustre/scratch119/humgen/projects/sc-eqtl-ibd/data/marker_gene_db-raw_data/database/celltypes/colon-immune/JamesKR-32066951/database.tsv' ]
            ]
        }




    }  
    bbknn{
        run_process = true
    }

    celltype_assignment{
        run_celltype_assignment=false
        run_azimuth=true
        run_keras=false
        run_celltypist=true
    }
    reduced_dims{
        vars_to_regress{
            value = ''   
        }
    }

}

process {

    withName: PLOT_DISTRIBUTIONS{
        containerOptions = "--containall --cleanenv --workdir /tmp -B /tmp"
    }

    withName: CELLEX_CLUSTER_MARKERS{
        maxForks=7
        memory = 300.GB
    }
    
    withName: GATHER_DATA{
        maxForks=7
        memory = 100.GB
    }
    withName: LISI{
        maxForks=7
        memory = 300.GB
    }
    withName: CLUSTER_VALIDATE_RESOLUTION_KERAS{
        memory = 300.GB
    }

    withName: UMAP_CALCULATE_AND_PLOT{
        memory = 300.GB
    }

    withName: SCCAF_ASSESS_CLUSTERING{
        memory = 300.GB
    }
    
}
