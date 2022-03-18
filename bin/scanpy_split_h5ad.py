#!/usr/bin/env python3

__date__ = '2020-03-13'
__version__ = '0.0.1'

import sys
import os
os.environ['NUMBA_CACHE_DIR']='/tmp'
os.environ['MPLCONFIGDIR']='/tmp'
from distutils.version import LooseVersion
import scipy
import scipy.io
import gzip
import pandas
import scanpy
import anndata
## split h5ad file by chromium channel

# for testing use
# infnam = "/lustre/scratch123/hgi/mdt1/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/franke_data/work/2a/ebdc5079ae949777263ce3b1aca510/BF61CE54F4603C9F-adata.h5ad"

def split_h5ad_by_batch(ad, oufnprfx, colnam_batch = 'batch', anndata_compression_opts=None, option='true'):
    oufn_list_fnam = '{}_files.txt'.format(oufnprfx)
    oufn_list = []
    batch_labels = pandas.Categorical(ad.obs[colnam_batch].apply(lambda a: a.split('__')[0])) # <class 'pandas.core.series.Series'>
    samples = {}
    count=0
    if (option == 'true'):
        print("Here we split the h5ad we just normalise for celltype assignmet")
        for bl in batch_labels.categories:
            print(bl)
            count+=1
            oufnam = '{0}_{1}.h5ad'.format(oufnprfx, bl)
            samples[count]={'experiment_id':bl,'h5ad_filepath':f"{os.getcwd()}/{oufnam}"}
            oufn_list.append(oufnam)
            adb = ad[batch_labels == bl,:]
            # strip unneccessary meta data - this seems to aid subsequent transformation to Seurat h5 format
            # assumes that cells failing qc already were stripped out
            # ad.obs = ad.obs[['convoluted_samplename', 'cell_passes_qc']]
            # ad.var = ad.var[['feature_types', 'genome', 'gene_symbols']]
            adb.obs = pandas.DataFrame(adb.obs.index, index = adb.obs.index, columns = ["cell_barcode"])

            try:
                vdf = adb.var[["feature_types", "genome"]]
            except:
                adb.var['genome']='GRCh38'
                vdf = adb.var[["feature_types", "genome"]]
            vdf.insert(1,"gene_ids", vdf.index)
            vdf.index = pandas.Index(adb.var['gene_symbols'].astype('str'))
            #ad.var = vdf.set_index("gene_symbols", drop = True, verify_integrity = False)
            adb.var = vdf

            del adb.uns
            adata = scanpy.AnnData(adb.X)
            adata.obs= adb.obs
            adata.var = adb.var
            adb = adata
            if anndata_compression_opts is None:
                adb.write(oufnam)
            else:
                adb.write(
                    oufnam,
                    compression='gzip',
                    compression_opts=anndata_compression_opts
                )
            oufn_list.append(oufnam)
    else:
        print("Here we dont split the h5ad we just normalise for celltype assignmet")
        bl ='full'
        count+=1
        oufnam = '{0}_{1}.h5ad'.format(oufnprfx, bl)
        samples[count]={'experiment_id':bl,'h5ad_filepath':f"{os.getcwd()}/{oufnam}"}
        adb = ad
        adb.obs = pandas.DataFrame(adb.obs.index, index = adb.obs.index, columns = ["cell_barcode"])
        try:
            vdf = adb.var[["feature_types", "genome"]]
        except:
            adb.var['genome']='GRCh38'
            vdf = adb.var[["feature_types", "genome"]]

        vdf.insert(1,"gene_ids", vdf.index)
        vdf.index = pandas.Index(adb.var['gene_symbols'].astype('str'))
        # vdf.index.name = None
        #ad.var = vdf.set_index("gene_symbols", drop = True, verify_integrity = False)
        adb.var = vdf
        del adb.layers
        adata = scanpy.AnnData(adb.X)
        adata.obs= adb.obs
        adata.var = adb.var
        adb = adata
        if anndata_compression_opts is None:
            adb.write(oufnam)
        else:
            adb.write(
                oufnam,
                compression='gzip',
                compression_opts=anndata_compression_opts)
        oufn_list.append(oufnam)

    Samples = pandas.DataFrame(samples).T
    Samples.to_csv('Samples.tsv',sep='\t',index=False)
    with open(oufn_list_fnam, 'w') as oufh:
        for fn in oufn_list:
            oufh.write(fn + '\n')
    return oufn_list_fnam

if __name__ == '__main__':
    nargs = len(sys.argv)
    # if nargs != 3:
    #     sys.exit("usage: %s ./<input_file *.h5ad> <output_file_prefix>".format(sys.argv[0]))
    infnam = sys.argv[1]
    oufnprfx = sys.argv[2]
    option = sys.argv[3]

    # print(infnam)
    ad = scanpy.read(infnam)
    split_h5ad_by_batch(ad, oufnprfx, colnam_batch = 'batch', anndata_compression_opts = None,option=option)
    sys.exit()
