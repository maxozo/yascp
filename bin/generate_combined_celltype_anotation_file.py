#!/usr/bin/env python

__date__ = '2021-11-15'
__version__ = '0.0.1'

import argparse
from distutils.version import LooseVersion
import os
os.environ['NUMBA_CACHE_DIR']='/tmp'
os.environ['MPLCONFIGDIR']='/tmp'
import pandas as pd
import scanpy


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Merges all the Celltype information in one file for Azimuth and Celltypist.
            """
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )

    parser.add_argument(
        '-az', '--all_azimuth_files',
        action='store',
        dest='all_azimuth_files',
        required=True,
        help='List of csv-delimited files of celltypes assigned by azimuth.'
    )

    parser.add_argument(
        '-ad', '--adata',
        action='store',
        dest='andata',
        required=True,
        help='Input adata to add labels to'
    )

    parser.add_argument(
        '-ct', '--all_celltypist_files',
        action='store',
        dest='all_celltypist_files',
        required=True,
        help='String of labels for each reduced_dims_tsv file. List should be\
            split by "::".'
    )

    options = parser.parse_args()
    azimuth_files = options.all_azimuth_files.split('::')
    
    Data_All=pd.DataFrame()
    for azimuth_file1 in azimuth_files:
        Data=pd.read_csv(azimuth_file1,compression='gzip',sep='\t')
        Data = Data.rename(columns={'predicted.celltype.l2':'Azimuth:predicted.celltype.l2','predicted.celltype.l2.score':'Azimuth:predicted.celltype.l2.score','mapping.score':'Azimuth:mapping.score'})
        Data_All = pd.concat([Data_All,Data])

    celltypist_files = options.all_celltypist_files.split('::')
    for celltypist_file1 in celltypist_files:
        if 'input.' in celltypist_file1:
            print('skip')
        else:
            Model = celltypist_file1.split('___')[1]
            Data=pd.read_csv(celltypist_file1,index_col=0)
            try:
                Data_All.loc[Data['predicted_labels'].index,f'Celltypist:{Model}']=Data['predicted_labels']
            except:
                Data_All[f'Celltypist:{Model}']=''
                Data_All.loc[Data['predicted_labels'].index,f'Celltypist:{Model}']=Data['predicted_labels']
    Data_All.to_csv('All_Celltype_Assignments.csv',sep='\t')

    adata = options.andata
    ad = scanpy.read(adata)
    for col in Data_All.columns:
        ad.obs[col]=Data_All[col]
    ad.write(
        'adata.h5ad',
        compression='gzip',
        compression_opts=9  # takes ages, but we want a small file for system
    )



if __name__ == '__main__':
    main()
