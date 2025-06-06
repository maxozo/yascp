#!/usr/bin/env python3

__date__ = '2021-01-15'
__version__ = '0.0.1'
# for help, run: python3 split_h5ad_per_donor.py --help
# command with example inputs arguments: python3 split_h5ad_per_donor.py --samplename ukbb_scrna9479582 --vireo_donor_ids_tsv /lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/nf_ci_scrna_deconv/results/vireo/vireo_ukbb_scrna9479582/donor_ids.tsv --filtered_matrix_h5 /lustre/scratch123/hgi/mdt1/projects/ukbb_scrna/pipelines/fetch_Submission_Data_Pilot_UKB/nextflow_ci/pipelines/../../results/iget_study_cellranger/5933/ukbb_scrna9479582/cellranger_ukbb_scrna9479582/filtered_feature_bc_matrix.h5

# import python libraries:
# on farm5, these libraries are installed in a conda environment: conda activate nextflow
import logging
import click
import sys
import argparse
import os
os.environ['NUMBA_CACHE_DIR']='/tmp'
os.environ['MPLCONFIGDIR']='/tmp'
import csv
import random
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import seaborn as sns
import plotnine as plt9
from plotnine.ggplot import save_as_pdf_pages

# CLI arguments:
@click.command()

# required arguments:
@click.option('--vireo_donor_ids_tsv', required=False, type=str,
              help='path to donor_ids.tsv file, which is an output file of Vireo')

@click.option('--filtered_matrix_h5', required=True, type=click.Path(exists=True),
              help='path to filtered cells h5 data file created by cellranger count, which was deconvoluted by Vireo')

@click.option('--samplename', required=True, type=str,
              help='sample name of cellranger experiment deconvoluted by Vireo')

# arguments that are optional because they have default value:
@click.option('-o','--output_dir', default='./vireo_deconv_out', show_default=True, type=str,
              help='output directory for script output plots and files')

@click.option('-m','--print_modules_version', default=False, show_default=True, type=bool,
              help='True or False: whether to print version of all python module to file.')

@click.option('-p','--plot_n_cells_per_vireo_donor', default=True, show_default=True,  type=bool,
              help='True or False: whether to plot number of cells per deconvoluted donor to pdf in --output_dir')

@click.option('-g','--input_h5_genome_version', default="GRCh38", show_default=True, type=str,
              help='True or False: whether to write donor level scanpy hdf5 objects to dir --output_dir')

@click.option('-s','--scrublet', default="None", show_default=True, type=str,
              help='the path to the scrublet multiplet input')

@click.option('-w','--write_donor_level_filtered_cells_h5', default=True, show_default=True, type=bool,
              help='True or False: whether to write donor level scanpy hdf5 objects to dir --output_dir')

@click.option('-c','--anndata_compression_level', default=6,
              type=click.IntRange(1, 9, clamp=True), show_default=True,
              help='Gzip compression level for scanpy write of AnnData hdf5 objects. Integer in range 1 to 9')

@click.option('-d','--plotnine_dpi', default=100,
              type=click.IntRange(1, 1000, clamp=True), show_default=True,
              help='DPI pdf plots resolution for plotnine plots. Integer in range 1 to 1000')


def split_h5ad_per_donor(vireo_donor_ids_tsv, filtered_matrix_h5, samplename,
                         output_dir, print_modules_version, plot_n_cells_per_vireo_donor,
                         input_h5_genome_version,
                         write_donor_level_filtered_cells_h5, plotnine_dpi,
                         anndata_compression_level,scrublet):
    """split_h5ad_donor main script"""
    logging.info('running split_h5ad_per_donor() function..')

    # Set seed for reproducibility
    seed_value = 0
    # 0. Set `PYTHONHASHSEED` environment variable at a fixed value
    # os.environ['PYTHONHASHSEED']=str(seed_value)
    # 1. Set `python` built-in pseudo-random generator at a fixed value
    random.seed(seed_value)
    # 2. Set `numpy` pseudo-random generator at a fixed value
    np.random.seed(seed_value)
    sns.set(style='whitegrid')
    # Set the default dpi
    plt9.options.dpi = plotnine_dpi

    if not os.path.exists(output_dir):
        logging.info('creating directory ' + output_dir)
        os.makedirs(output_dir)

    # print modules version to output file.
    if print_modules_version:
        with open(output_dir + '/module_versions.txt', 'w') as f:
            for name, module in sorted(sys.modules.items()):
                if hasattr(module, '__version__'):
                    f.write(str(name) + '=' + str(module.__version__) + '\n')

    # fix genome key in cellbender h5 output
    # import tables
    
    orig_h5 = filtered_matrix_h5
    # the next bits are for reading in h5ad or h5 files instead of mtx
    # fixed_h5 = 'fixed_genome.h5'
    # tables.copy_file(orig_h5, fixed_h5, overwrite = True)
    # with tables.open_file(fixed_h5, "r+") as f:
    #     n = f.get_node("/matrix/features")
    #     n_genes = f.get_node("/matrix/shape")[0]
    #     if "genome" not in n:
    #         f.create_array(n, "genome", np.repeat(input_h5_genome_version, n_genes))
    # # read-in cellranger 10x data produced by 'cellranger count':
    # logging.info('fixed orig_h5 into fixed_h5')
    adata = sc.read_10x_mtx(orig_h5) #, genome='background_removed')
    # adata = sc.read_10x_h5(filtered_matrix_h5) #, genome='background_removed')

    logging.info(adata.var)
    logging.info(adata.obsm)
    logging.info("n cells len(adata.obs): " + str(len(adata.obs)))
    # logging.info('exiting..'); sys.exit()

    adata.var['gene_symbols'] = adata.var.index
    adata.var.index = adata.var['gene_ids'].values
    del adata.var['gene_ids']
    method = 'Vireo'
        
    # also read-in the cell deconvolution annotation produced by Vireo:
    if vireo_donor_ids_tsv !='None':
        vireo_anno_deconv_cells = pd.read_csv(vireo_donor_ids_tsv, sep='\t',
                                            index_col='cell')

        # calculate number of cells per deconvoluted donor:
        cells_per_donor_count = vireo_anno_deconv_cells[['donor_id']].value_counts().to_frame('n_cells')
        cells_per_donor_count.reset_index(level=cells_per_donor_count.index.names, inplace=True)
        logging.info('cells_per_donor_count:')
        logging.info(cells_per_donor_count)
        logging.info('sum(cells_per_donor_count.n_cells: ' + str(sum(cells_per_donor_count.n_cells)))

        # check that `adata` and `vireo_anno_deconv_cells` indexes DO match, as expected:
        logging.info('len(vireo_anno_deconv_cells.index): ' + str(len(vireo_anno_deconv_cells.index)))
        for cell_in_vireo_index in vireo_anno_deconv_cells.index:
            if cell_in_vireo_index not in adata.obs_names:
                logging.info('warning: cell index ' + cell_in_vireo_index + ' is in vireo_anno_deconv_cells but not in adata')

        # add Vireo annotation to adata
        adata.obs['convoluted_samplename'] = samplename
        for new_cell_annotation in ['donor_id','prob_max','prob_doublet','n_vars','best_singlet','best_doublet']:
            if new_cell_annotation in vireo_anno_deconv_cells.columns:
                logging.info('adding vireo annotation ' + new_cell_annotation + ' to AnnData object.')
                adata.obs[new_cell_annotation] = vireo_anno_deconv_cells[new_cell_annotation]
            else:
                adata.obs[new_cell_annotation] = 'nan'
                logging.info('warning: column ' + new_cell_annotation + ' is not in input Vireo annotation tsv.')
    else:
        logging.info('Samples are not deconvoluted')
        adata.obs['convoluted_samplename'] = samplename
        method = 'After cellbender'
        donor_cell_nr = len(adata.obs)
        cells_per_donor_count_dic = {}
        cells_per_donor_count_dic[1]={"donor_id":'donor','n_cells':donor_cell_nr}
        adata.obs['donor_id']='donor'
        # cells_per_donor_count_dic[2]={'donor_id':'doblets','n_cells':doublets_nr}
        cells_per_donor_count = pd.DataFrame(cells_per_donor_count_dic).T
    
    # plot n cells per deconvoluted Vireo donor:
    if plot_n_cells_per_vireo_donor:
        gplt = plt9.ggplot(cells_per_donor_count, plt9.aes(
            x='donor_id',
            y='n_cells',
            fill='donor_id'
        ))
        gplt = gplt + plt9.theme_bw() + plt9.theme(legend_position='none',
                                                   axis_text_x=plt9.element_text(colour="black", angle=45),
                                                   axis_text_y=plt9.element_text(colour="black"))
        gplt = gplt + plt9.geom_bar(stat='identity', position='dodge')
        gplt = gplt + plt9.geom_text(plt9.aes(label='n_cells'))
        gplt = gplt + plt9.labels.ggtitle(f'{method} \nnumber of cells per deconvoluted donor\nsample: ' + samplename)
        gplt = gplt + plt9.labels.xlab('deconvoluted donor')
        gplt = gplt + plt9.labels.ylab(f'Number of cells assigned by {method}')

        # save plot(s) as pdf:
        plots = [gplt]
        save_as_pdf_pages(plots, filename=output_dir + '/Vireo_plots.pdf')

    # write AnnData with Vireo cell annotation in .obs
    output_file = output_dir + '/vireo_annot.' + samplename
    logging.info('Write h5ad AnnData with Vireo cell annotation in .obs to ' + output_file)
    adata.write('{}.h5ad'.format(output_file), compression='gzip', compression_opts= anndata_compression_level)

    if write_donor_level_filtered_cells_h5:

        if not os.path.exists(output_dir + '/donor_level_anndata'):
            logging.info('creating directory ' + output_dir + '/donor_level_anndata')
            os.makedirs(output_dir + '/donor_level_anndata')

        adata_donors = []
        adata_nr_cells = {}
        count=0
        for donor_id in adata.obs['donor_id'].unique():
            donor_id = str(donor_id)
            logging.info('filtering cells of AnnData to donor ' + donor_id)
            adata_donor = adata[adata.obs['donor_id'] == donor_id, :]
            
            logging.info("n cells len(adata_donor.obs) for " + donor_id  + ': ' + str(len(adata_donor.obs)) + '/' + str(len(adata.obs)))
            if len(adata_donor.obs) > 0:
                
                count+=1
                logging.info("more than 0 cells for donor")
                adata_nr_cells[count]={'experiment_id':f'{samplename}__{donor_id}','n_cells':len(adata_donor.obs)}
                adata_donors.append((donor_id, adata_donor))
                output_file = output_dir + '/donor_level_anndata/' + donor_id + '.' + samplename
                pd.DataFrame(adata_donor.obs.index,columns=['barcodes']).to_csv(f'{output_file}.barcodes.tsv',index=False,header=False)
                logging.info('Write h5ad donor AnnData to ' + output_file)
                adata_donor.write('{}.h5ad'.format(output_file), compression='gzip', compression_opts= anndata_compression_level)
            else:
                logging.info("0 cells for donor, therefore no writing donor-specific h5ad.")
        adata_nr_cells_df = pd.DataFrame(adata_nr_cells).T
        adata_nr_cells_df.to_csv('exp__donor_n_cells.tsv',sep='\t', index=False,header=False)

if __name__ == '__main__':
    # set logging level and handler:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s",
                        handlers=[logging.StreamHandler()]) # logging.FileHandler("debug.log"),
    split_h5ad_per_donor()

