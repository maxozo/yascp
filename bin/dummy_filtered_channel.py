#!/usr/bin/env python3

__date__ = '2021-11-04'
__version__ = '0.0.1'

import argparse
import os
os.environ['NUMBA_CACHE_DIR']='/tmp'
os.environ['MPLCONFIGDIR']='/tmp'
import scanpy as sc
import pandas as pd


def main(ad_in=False,id_in='sanger_sample_id'):
    """
        prepearinf file like: 
            "experiment_id" "filter_type"   "n_cells_left_in_adata"
            "franke_Pilot_3_lane_1__s1"     "before_filters"        604
            "franke_Pilot_3_lane_1__s2"     "before_filters"        774
            "franke_Pilot_3_lane_1__s3"     "after_filters" 898
            "franke_Pilot_3_lane_1__s4"     "after_filters" 615
    """

    print("Creating dummy filtered file")
    Adata_Object = sc.read_h5ad(filename=ad_in)
    data = {}
    count=1
    for id in set(Adata_Object.obs[id_in]):
        count+=1
        ncells = len(Adata_Object[Adata_Object.obs[id_in]==id])
        data[count]={'experiment_id':id,'filter_type':"before_filters",'n_cells_left_in_adata':ncells}
        data[count]={'experiment_id':id,'filter_type':"after_filters",'n_cells_left_in_adata':ncells}
    d2 = pd.DataFrame(data).T
    d2.to_csv('filtered_cell_dummy.tsv',sep='\t',index=False)
    print("Creation finished")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""
            Read anndata object. Runs BBKNN. Assumes PCS have already been
            calculated.
            """
    )
    parser.add_argument(
        '-h5', '--h5_anndata',
        action='store',
        dest='h5',
        required=True,
        help='H5 AnnData file.'
    )

    parser.add_argument(
        '-id', '--id',
        action='store',
        dest='id',
        required=True,
        help='H5 AnnData id collum.'
    )

    options = parser.parse_args()
    ad_in = options.h5
    id_in = options.id
    main(ad_in=ad_in,id_in=id_in)
