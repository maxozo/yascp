#!/usr/bin/env python

__date__ = '2020-07-09'
__version__ = '0.0.1'

import argparse
import os
os.environ['NUMBA_CACHE_DIR']='/tmp'
os.environ['MPLCONFIGDIR']='/tmp'
from distutils.version import LooseVersion
from typing import Dict
import tables
import scipy
import scipy.io
import scipy.sparse
import gzip
import pandas as pd
import numpy as np
import anndata
filter_0_count_cells=True

def dict_from_h5(file: str) -> Dict[str, np.ndarray]:
    """Read in everything from an h5 file and put into a dictionary."""
    d = {}
    with tables.open_file(file) as f:
        # read in everything
        for array in f.walk_nodes("/", "Array"):
            d[array.name] = array.read()
    return d


def anndata_from_h5(file: str,analyzed_barcodes_only: bool = True) -> 'anndata.AnnData':
    """Load an output h5 file into an AnnData object for downstream work.

    Args:
        file: The h5 file
        analyzed_barcodes_only: False to load all barcodes, so that the size of
            the AnnData object will match the size of the input raw count
            matrix. True to load a limited set of barcodes: only those
            analyzed by the algorithm. This allows relevant latent
            variables to be loaded properly into adata.obs and adata.obsm,
            rather than adata.uns.

    Returns:
        adata: The anndata object, populated with inferred latent variables
            and metadata.
    """
    d = dict_from_h5(file)
    X = scipy.sparse.csc_matrix(
        (d.pop('data'), d.pop('indices'), d.pop('indptr')),
        shape=d.pop('shape')
    ).transpose().tocsr()

    if analyzed_barcodes_only:
        if 'barcodes_analyzed_inds' in d.keys():
            X = X[d['barcodes_analyzed_inds'], :]
            d['barcodes'] = d['barcodes'][d['barcodes_analyzed_inds']]
        elif 'barcode_indices_for_latents' in d.keys():
            X = X[d['barcode_indices_for_latents'], :]
            d['barcodes'] = d['barcodes'][d['barcode_indices_for_latents']]
        else:
            print(
                'Warning: analyzed_barcodes_only=True, but the key ',
                '"barcodes_analyzed_inds" or "barcode_indices_for_latents" ',
                'is missing from the h5 file. ',
                'Will output all barcodes, and proceed as if ',
                'analyzed_barcodes_only=False'
            )

    print(d.keys())

    # Construct the count matrix.
    if 'gene_names' in d.keys():
        gene_symbols = d.pop('gene_names').astype(str)
    else:
        gene_symbols = d.pop('name').astype(str)
    adata = anndata.AnnData(
        X=X,
        obs={'barcode': d.pop('barcodes').astype(str)},
        var={
            'gene_ids': d.pop('id').astype(str),
            'gene_symbols': gene_symbols,
            'feature_type':d.pop('feature_type').astype(str),
        }
    )
    # adata = adata[:,adata.var.query('feature_type=="Gene Expression"').index]
    adata.obs.set_index('barcode', inplace=True)
    adata.var.set_index('gene_ids', inplace=True)
    
    
    
    # Add other information to the adata object in the appropriate slot.
    for key, value in d.items():
        try:
            value = np.asarray(value)
            if len(value.shape) == 0:
                adata.uns[key] = value
            elif value.shape[0] == X.shape[0]:
                if (len(value.shape) < 2) or (value.shape[1] < 2):
                    adata.obs[key] = value
                else:
                    adata.obsm[key] = value
            elif value.shape[0] == X.shape[1]:
                if value.dtype.name.startswith('bytes'):
                    adata.var[key] = value.astype(str)
                else:
                    adata.var[key] = value
            else:
                adata.uns[key] = value
        except Exception:
            print(
                'Unable to load data into AnnData: ', key, value, type(value)
            )

    if analyzed_barcodes_only:
        cols = adata.obs.columns[
            adata.obs.columns.str.startswith('barcodes_analyzed')
            | adata.obs.columns.str.startswith('barcode_indices')
        ]
        for col in cols:
            try:
                del adata.obs[col]
            except Exception:
                pass
            
    if adata.obs_names.duplicated().any():
        print("Warning: duplicated barcodes found — making obs_names unique.")
        adata.obs_names_make_unique()
    return adata



def cellbender_to_tenxmatrix(adata,out_file='',out_dir='tenx_from_adata',verbose=True):
    """Write 10x like data from 10x H5.

    Parameters
    ----------
    adata : pandas.DataFrame
        Description of parameter `adata`.
    output_file : string
        Description of parameter `output_file`.
    out_dir : string
        Description of parameter `out_dir`.
    verbose : boolean
        Description of parameter `verbose`.

    Returns
    -------
    execution_code : int
    """
    # Make the output directory if it does not exist.
    if out_dir == '':
        out_dir = os.getcwd()
    else:
        os.makedirs(out_dir, exist_ok=True)

    import pandas as pd

    # Suppose adata is your AnnData object.
    # Define a function to check if a string starts with 'ENSG'
    def looks_like_ensembl(id_str):
        # A simple heuristic: Ensembl gene IDs typically start with 'ENSG'
        return isinstance(id_str, str) and id_str.startswith("ENSG")

    # Get the gene symbols column as a Pandas Series
    gene_symbols_series = adata.var['gene_symbols'].astype(str)

    # Calculate the proportion of values that look like an Ensembl ID in the gene_symbols column.
    prop_in_gene_symbols = gene_symbols_series.str.startswith("ENSG").mean()

    # Calculate the proportion in the index (adata.var_names)
    prop_in_index = sum(looks_like_ensembl(gene) for gene in adata.var_names) / len(adata.var_names)

    print("Proportion of ENSG pattern in gene_symbols column: {:.2f}".format(prop_in_gene_symbols))
    print("Proportion of ENSG pattern in index: {:.2f}".format(prop_in_index))

    # -----------------------------------------
    # Determine gene_ids and gene_symbols
    # -----------------------------------------
    threshold = 0.5

    if 'gene_symbols' in adata.var.columns:
        gene_symbols_series = adata.var['gene_symbols'].astype(str)
    else:
        gene_symbols_series = pd.Index(adata.var_names).astype(str)

    def looks_like_ensembl(x):
        return isinstance(x, str) and x.startswith("ENSG")

    prop_in_index = pd.Series(adata.var_names).map(looks_like_ensembl).mean()
    prop_in_gene_symbols = gene_symbols_series.map(looks_like_ensembl).mean()

    if prop_in_gene_symbols >= threshold and prop_in_index < threshold:
        print("Using index as gene symbols and gene_symbols column as Ensembl IDs.")
        adata.var['gene_ids'] = gene_symbols_series
        # index stays as is (gene symbols)

    elif prop_in_index >= threshold and prop_in_gene_symbols < threshold:
        print("Using gene_symbols column as gene symbols and index as Ensembl IDs.")
        adata.var['gene_ids'] = adata.var_names
        adata.var_names = gene_symbols_series

    elif prop_in_index >= threshold and prop_in_gene_symbols >= threshold:
        print("Warning: Both index and gene_symbols look like Ensembl IDs.")
        adata.var['gene_ids'] = adata.var_names  # fallback
    else:
        print("Warning: Neither look clearly like Ensembl IDs.")
        adata.var['gene_ids'] = adata.var_names  # fallback

    out_f = os.path.join(out_dir, f"{out_file}barcodes.tsv.gz")
    if verbose:
        print(f'Writing {out_f}')
    pd.DataFrame(adata.obs.index).to_csv(
        out_f,
        sep='\t',
        index=False,
        header=False,
        compression='gzip'
    )

    # -----------------------------------------
    # Create features.tsv.gz content
    # -----------------------------------------
    out_f = os.path.join(out_dir, f"{out_file}features.tsv.gz")
    if verbose:
        print(f'Writing {out_f}')

    try:
        df_features = pd.DataFrame(
            data=[
                adata.var['gene_ids'].values,
                adata.var_names.values,
                adata.var.feature_type.values
            ]
        ).T
    except AttributeError:
        df_features = pd.DataFrame(
            data=[
                adata.var['gene_ids'].values,
                adata.var_names.values,
                adata.var.feature_types.values
            ]
        ).T


    df_features.to_csv(
        out_f,
        sep='\t',
        compression='gzip',
        index=False,
        header=False
    )

    # Save the count-adjusted matrix
    out_mtx = adata.X.transpose()
    out_mtx = out_mtx.astype(np.int32)
    if not isinstance(out_mtx, scipy.sparse.csr.csr_matrix):
        out_mtx = scipy.sparse.csr_matrix(out_mtx)
    out_f = os.path.join(
        out_dir,
        '{}matrix.mtx.gz'.format(out_file)
    )
    if verbose:
        print('Writing {}'.format(out_f))
        # print(out_mtx.dtype)  # it looks like even counts stored as float
    with gzip.open(out_f, 'wb', compresslevel=9) as f:
        scipy.io.mmwrite(
            f,
            out_mtx,
            comment='metadata_json: {"format_version": 2}'
            # field='real'  # can be integer if counts otherwise real
        )

    if verbose:
        print(
            'Done.'
        )

    return 0


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Read an H5 10x-formatted object and write matrix files similar to
            10x output. This script requires pandas >1.0.1
            """
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )

    parser.add_argument(
        '--nf_outdir_tag',
        action='store',
        dest='nf_outdir_tag',
        required=True,
        help='Output directory tag used in nextflow.'
    )

    parser.add_argument(
        '--cb_outfile_tag',
        action='store',
        dest='cb_outfile_tag',
        required=True,
        help='Output file tag used in CellBender.'
    )

    parser.add_argument(
        '--experiment_id',
        action='store',
        dest='experiment_id',
        required=True,
        help='Sample experiment id.'
    )

    parser.add_argument(
        '-f', '--fpr',
        action='store',
        dest='fpr',
        required=True,
        help='False postive rate tag used in CellBender call.'
    )

    parser.add_argument(
        '--cb_params',
        action='store',
        dest='cb_params',
        required=True,
        help='Cellbender parameters.'
    )

    # parser.add_argument(
    #     '-of', '--output_file',
    #     action='store',
    #     dest='output_file',
    #     default='cellbender_results',
    #     help='Basename of output files.\
    #         (default: %(default)s)'
    # )

    options = parser.parse_args()

    # Get a list of all of the output files
    cellbender_out_files = {}
    fprs = options.fpr.split(' ')
    if len(fprs) == 1:  # then no FPR in the output file
        file = '{}_filtered.h5'.format(options.cb_outfile_tag)
        cellbender_out_files[fprs[0]] = file
    else:
        for fpr in fprs:
            file = '{}_FPR_{}_filtered.h5'.format(
                options.cb_outfile_tag,
                fpr.replace('.', 'pt')
            )
            cellbender_out_files[fpr] = file

    res_list = []
    
    for fpr, fil in cellbender_out_files.items():
        # Check to make sure each file exists
        if not os.path.isfile(fil):
            # File may not exist because of prior edit to replace . with pt
            fil = fil.replace('pt', '.')
            if not os.path.isfile(fil):
                raise Exception('File does not exist:\t{}'.format(fil))

        # Move the unfiltered files
        try:
            os.rename(
                fil.replace('_filtered', ''),
                fil.replace('_filtered', '_unfiltered')
            )
        except:
            print('Already_moved')

        local_outdir = '{}-FPR_{}-filtered_10x_mtx'.format(
            options.cb_outfile_tag,
            fpr.replace('.', 'pt')
        )
        
        local_outdir_raw = '{}-FPR_{}-unfiltered_10x_mtx'.format(
            options.cb_outfile_tag,
            fpr.replace('.', 'pt')
        )        
        
        print(local_outdir, fil, local_outdir)
        res_list.append({
            'fpr': fpr,
            'cb_params': options.cb_params,
            'experiment_id': options.experiment_id,
            'data_path_10x_format': local_outdir
        })
        # Load the AnnData file
        # adata = sc.read_10x_h5(filename=options.h5, genome=options.g)
        # Load remove-background output data.
        # We need to do this because of bug here:
        # https://github.com/broadinstitute/CellBender/issues/57

        adata_raw = anndata_from_h5(
            fil.replace('_filtered', '_unfiltered'),
            analyzed_barcodes_only=False
        )

        if (adata_raw.X.data < 0).any():
            adata_raw.X.data[adata_raw.X.data < 0] = 0
            print("Cellbender produced negative values which we replace with 0s.")
            with open("Warnings.log", "a") as log_file:  # "a" mode appends to the file
                log_file.write(f"Cellbender {fpr} unfiltered produced negative values which we replace with 0s.\n")            
        # Run the conversion function.
        _ = cellbender_to_tenxmatrix(
            adata_raw,
            out_file='',
            out_dir='{}-FPR_{}-unfiltered_10x_mtx'.format(
                options.cb_outfile_tag,
                fpr.replace('.', 'pt')
            )
        )
        del adata_raw
          
        try:
            adata = anndata_from_h5(
                fil,
                analyzed_barcodes_only=True
            )
        except:
            # issue is fixed in the new cellbender version so this can be loaded differently
            adata = anndata_from_h5(
                fil,
                analyzed_barcodes_only=False
            )
            
        if (adata.X.data < 0).any():
            adata.X.data[adata.X.data < 0] = 0
            print("Cellbender produced negative values which we replace with 0s.")
            with open("Warnings.log", "a") as log_file:  # "a" mode appends to the file
                log_file.write(f"Cellbender {fpr} filtered produced negative values which we replace with 0s.\n")
            
            
        _ = cellbender_to_tenxmatrix(
            adata,
            out_file='',
            out_dir='{}-FPR_{}-filtered_10x_mtx'.format(
                options.cb_outfile_tag,
                fpr.replace('.', 'pt')
            )
        )



    # Now write the output file
    df = pd.DataFrame(res_list)
    nf_outdir_tag = options.nf_outdir_tag.rstrip('/')
    if nf_outdir_tag != '':
        df['data_path_10x_format'] = nf_outdir_tag + '/' + \
            df['data_path_10x_format']

    df.to_csv(
        '{}-filtered_10x_mtx-file_list.tsv'.format(options.cb_outfile_tag),
        sep='\t',
        index=False,
        header=True
    )


if __name__ == '__main__':
    main()
