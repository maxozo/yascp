#!/usr/bin/env python3

__date__ = '2021-18-11'
__version__ = '0.0.1'


import os
os.environ['NUMBA_CACHE_DIR']='/tmp'
os.environ['MPLCONFIGDIR']='/tmp'
import sys
import pandas as pd
import scanpy
import argparse

def main(andata,anndata_compression_level,metadata):
    print('lets add the extra metadata to andata')
    print(andata)
    adata = scanpy.read_h5ad(filename=andata)
    
    metadata_data = pd.read_csv(metadata,sep='\t',index_col='experiment_id')
    for col1 in metadata_data.columns:
        print(col1)
        try:
            adata.obs[col1]=''
            for experiment_id in metadata_data.index:
                
                metadata_val = metadata_data.loc[experiment_id,col1]
                adata.obs.loc[adata.obs['convoluted_samplename']==experiment_id,col1]=metadata_val
                print(f"{experiment_id} val : {metadata_val}")           
        except:
            sys.exit(f'Metadata column {col1} is already in andata, please rename or remove this from metadata file!')
        
                # t = adata.obs

    adata.write(
    'andata_with_metadata.h5ad',
    compression='gzip',
    compression_opts=int(anndata_compression_level))


if __name__ == '__main__':
    
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
        '-ad', '--andata',
        action='store',
        dest='andata',
        required=True,
        help='Input adata to add labels to'
    )

    parser.add_argument(
        '-c', '--anndata_compression_level',
        action='store',
        dest='anndata_compression_level',
        required=True,
        help='Gzip compression level for scanpy write of AnnData hdf5 objects. Integer in range 1 to 9'
    )

    parser.add_argument(
        '-m', '--metadata',
        action='store',
        dest='metadata',
        required=True,
        help='Tsv file containing metadata to be added to andata object.'
    )


    options = parser.parse_args()
    anndata_compression_level=options.anndata_compression_level
    andata=options.andata
    metadata=options.metadata
    main(andata,anndata_compression_level,metadata)



