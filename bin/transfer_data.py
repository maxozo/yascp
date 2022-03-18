#!/usr/bin/env python

__date__ = '2021-07-28'
__version__ = '0.0.1'

import os
from shutil import copyfile,copytree
from os import listdir
import glob
import argparse
import pandas as pd



def main_data_colection(pipeline='',name='',directory='',input_table=None,cb_res=None):
    
    name_dir=f"Summary_plots/{name}"
    try:
        os.makedirs(f'{name_dir}')
    except:
        print("dir exists")


    #Cellbender
    print('prepearing Cellbender folder')
    try:
        os.mkdir(f'{name_dir}/Cellbender')
    except:
        print('dire exists')
    
    if (cellbender)=='cellranger':
        # here we do not use cellbender and go with default cellranger
        df_cellbender = None
    elif (cellbender)=='cellbender':
        # here we have run the cellbender as par of pipeline. 
        file_path = glob.glob(f'{results_dir}/nf-preprocessing/cellbender/qc_cluster_input_files/*{cb_res}*')[0]
        df_cellbender = pd.read_table(file_path, index_col = 'experiment_id')
    else:
        # this is existing_cellbender, hence using this input
        df_cellbender = pd.read_table(f'{cellbender}', index_col = 'experiment_id')
    
    
    # folder1 = f'{directory}/nf-preprocessing/cellbender'
    # folder2 = f'{directory}/cellbender_vs_cellranger'
    
    resolution2=cb_res.replace('pt','.')
    if (df_cellbender is not None):
        print('t')
        for folder in df_cellbender.index:
            print(folder)
            dir1 = f"{df_cellbender.loc[folder,'data_path_10x_format']}/.."
            dir = f"{df_cellbender.loc[folder,'data_path_10x_format']}/../.."
            if os.path.isdir(dir1):
                copyfile(f'{dir1}/plots/cellbender_results-cellbender_FPR_{cb_res}_filtered-ambient_signature-scatter_genenames.png', f'{name_dir}/Cellbender/{folder}_ambient_signature-scatter_genenames.png')
                copyfile(f'{dir1}/plots/cellbender_results-cellbender_FPR_{cb_res}_filtered-abs_count_difference-boxplot.png', f'{name_dir}/Cellbender/{folder}_ount_difference-boxplot.png')
                copyfile(f'{dir1}/plots/cellbender.pdf', f'{name_dir}/Cellbender/cellbender_{folder}.pdf')
                copyfile(f'{dir}/compare_cellranger/fpr_{resolution2}/boxplots_cellranger_vs_cellbender.png', f'{name_dir}/Cellbender/{folder}_boxplots_cellranger_vs_cellbender.png')
                copyfile(f'{dir}/compare_cellranger/fpr_{resolution2}/barcode_vs_total_counts.png', f'{name_dir}/Cellbender/{folder}_barcode_vs_total_counts.png')
                copyfile(f'{dir}/compare_cellranger/fpr_{resolution2}/boxplot_topgenes_cellranger_vs_cellbender.png', f'{name_dir}/Cellbender/{folder}_boxplot_topgenes_cellranger_vs_cellbender.png')

    # Fetch Gather
    df_raw = pd.read_table(input_table, index_col = 'experiment_id')
    # 'Note that the names for the future projects may be different - have to be handled on the Nextflow modules'
    print('prepearing fetch folder')
    try:
        os.mkdir(f'{name_dir}/Fetch Pipeline')
    except:
        print('exists')
    metadata_table = pd.DataFrame()
    for folder in df_raw.index:
        dir3 = f"{df_raw.loc[folder, 'data_path_10x_format']}"
        copyfile(f'{dir3}/web_summary.html', f'{name_dir}/Fetch Pipeline/html_{folder}.html')
        metadata = pd.read_csv(f'{dir3}/metrics_summary.csv',sep=',',index_col=False)
        metadata['Sample_id']=folder
        metadata.set_index('Sample_id',drop=True,inplace=True)
        metadata_table=pd.concat([metadata_table,metadata])
    metadata_table.to_csv(f'{name_dir}/Fetch Pipeline/Submission_Data_Pilot_UKB.file_metadata.tsv')

    #if (pipeline=='Deconvolution'):
    folder1 = f'{directory}/deconvolution/split_donor_h5ad'
    if os.path.isdir(folder1):
        print('prepearing Deconvolution folder')
        try:
            os.mkdir(f'{name_dir}/Deconvolution')
        except:
            print('dire exists')
        Folders = listdir(folder1)
        for folder in Folders:
            copyfile(f'{folder1}/{folder}/Vireo_plots.pdf', f'{name_dir}/Deconvolution/Vireo_plots_{folder}.pdf')
        
    folder1 = f'{directory}/celltype/azimuth'
    if os.path.isdir(folder1):
        print('prepearing Deconvolution folder')
        try:
            os.mkdir(f'{name_dir}/Cell-type assignment')
        except:
            print('dire exists')
        Folders = listdir(folder1)
        for file in Folders:
            copyfile(f'{folder1}/{file}', f'{name_dir}/Cell-type assignment/{file}')

    folder1 = f'{directory}/plots'
    if os.path.isdir(folder1):
        try:
            os.mkdir(f'{name_dir}/QC metrics')
        except:
            print('dire exists')
        copyfile(f'{folder1}/adata-cell_desity.png', f'{name_dir}/QC metrics/adata-cell_desity.png')
        copyfile(f'{folder1}/adata-cell_filtered_per_experiment-n_cells_before_after.png', f'{name_dir}/QC metrics/adata-cell_filtered_per_experiment-n_cells_before_after.png')
        copyfile(f'{folder1}/scatterplot-sex_sample_swap_check.png', f'{name_dir}/QC metrics/scatterplot-sex_sample_swap_check.png')
        copyfile(f'{folder1}/adata-outlier_cells.png', f'{name_dir}/QC metrics/adata-outlier_cells.png')
        fil1 = glob.glob(f'{folder1}/plot_ecdf-x_log10*total_counts*')[0]
        copyfile(fil1, f'{name_dir}/QC metrics/plot_ecdf-x_log10.var=total_counts.color=experiment_id-adata.png')
        

    folder1 = f'{directory}/clustering'
    if os.path.isdir(folder1):
        try:
            os.mkdir(f'{name_dir}/Clustering')
        except:
            print('dire exists')
        Harmony_UMAPS = glob.glob(f'{folder1}/*/*harmony*/*/plots/umap*')
        for umap1 in Harmony_UMAPS:
            
            name = umap1.split('/')[-1]
            resolution = umap1.split('/')[-3].split("resolution=")[1]
            
            copyfile(umap1, f'{name_dir}/Clustering/res={resolution}_Harmony_{name}')
    
    folder1 = f'{directory}/UMAPs'
    if os.path.isdir(folder1):
        try:
            os.mkdir(f'{name_dir}/Clustering')
        except:
            print('dire exists')
        Coloured_UMAPS = glob.glob(f'{folder1}/*')
        for umap1 in Coloured_UMAPS:
            name = umap1.split('/')[-1]
            copyfile(umap1, f'{name_dir}/Clustering/{name}')                    
    
    
    folder1 = f'{directory}/handover/minimal_dataset_summary'
    if os.path.isdir(folder1):
        try:
            copytree(folder1, f'{name_dir}/Summary')
        except:
            print('dire exists')
        
        

    # if secret params exist then transfer data to web. 
    # rsync -vr Summary_plots ubuntu@1xxxxx:/volume/scRNA_test_app/scrna_static_and_media_files/media

if __name__ == "__main__":
    # This code will gather the most important graphs and place them in a summary folder.
    parser = argparse.ArgumentParser(
        description="""
            Collect important graphs for web visualisations
            """
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )

    parser.add_argument(
        '-rd', '--results_dir',
        action='store',
        dest='results_dir',
        required=True,
        help='The nextflow results directory.'
    )

    parser.add_argument(
        '-cb_res', '--cb_res',
        action='store',
        dest='cb_res',
        required=True,
        help='The cellbender resolution used in downstream analysis'
    )

    parser.add_argument(
        '-it', '--input_table',
        action='store',
        dest='input_table',
        required=True,
        help='Input table with paths to 10xRuns'
    )

    parser.add_argument(
        '-cb', '--cellbender',
        action='store',
        dest='cellbender',
        required=True,
        help='The path to cellbender runs, since this may be external.'
    )
    options = parser.parse_args()
    
    t=os.getcwd()
    results_dir = options.results_dir
    name = t.split('/')[-4]
    
    input_table=options.input_table
    cb_res=options.cb_res
    cellbender = options.cellbender
    main_data_colection(name=f"{name}",directory=results_dir,input_table=input_table,cb_res=cb_res)
