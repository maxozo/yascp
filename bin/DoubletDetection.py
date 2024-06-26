#!/usr/bin/env python
import numpy as np
import doubletdetection
import tarfile
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt
import os
import argparse
import sys
import pandas as pd
import scanpy as sc

# Load read10x function from mods directory

# mods_path = "/opt/WG1-pipeline-QC/Demultiplexing/mods"
# sys.path.append(mods_path)
# import read10x


parser = argparse.ArgumentParser(
    description="wrapper for DoubletDetection for doublet detection from transcriptomic data.")
# parser.add_argument("-m", "--counts_matrix", required = True, help = "cell ranger counts matrix.mtx")
# parser.add_argument("-b", "--barcodes", required = True, help = "cell ranger barcodes.tsv or barcodes.tsv.gz")
parser.add_argument("-o", "--outdir", required = False, default = os.getcwd(), help = "The output directory; default is current working directory")
parser.add_argument("-i", "--n_iterations", required = False, default = 50, type = int, help = "Number of iterations to use; default is 50")
parser.add_argument("-p", "--phenograph", required = False, default = False, help = "Whether to use phenograph (True) or not (False); default is False")
parser.add_argument("-s", "--standard_scaling", required = False, default = False, help = "Whether to use standard scaling of normalized count matrix prior to PCA (True) or not (False); default is True")
parser.add_argument("-t", "--p_thresh", required = False, default = 1e-16, type = float, help = "P-value threshold for doublet calling; default is 1e-16")
parser.add_argument("-v", "--voter_thresh", required = False, default = 0.5, type = float, help = "Voter threshold for doublet calling; default is 0.5")
parser.add_argument(
    '-txd', '--tenxdata_dir',
    action='store',
    dest='txd',
    default='',
    required=True,
    help='Path to directory with data in 10x matrix format.'
)
args = parser.parse_args()
# adata = sc.read_10x_h5(args.txd)
adata = sc.read_10x_mtx(
    path=args.txd,
    # var_names='gene_symbols',
    var_names='gene_ids',
    make_unique=False
)

print(args.p_thresh)
print(args.voter_thresh)
print(args.phenograph)
if args.phenograph == 'True':
    pheno = True
elif args.phenograph == 'False':
    pheno = False
else:
    pheno = args.phenograph
print(pheno)

if args.standard_scaling == 'True':
    standard_scaling = True
elif args.standard_scaling == 'False':
    standard_scaling = False
else:
    standard_scaling = args.standard_scaling
print(standard_scaling)


### Read in data ###
raw_counts =  adata.X
barcodes_df = adata.obs.index

print('Counts matrix shape: {} rows, {} columns'.format(raw_counts.shape[0], raw_counts.shape[1]))

# Remove columns with all 0s
zero_genes = (np.sum(raw_counts, axis=0) == 0).A.ravel()
raw_counts = raw_counts[:, ~zero_genes]
print('Counts matrix shape after removing unexpressed genes: {} rows, {} columns'.format(raw_counts.shape[0], raw_counts.shape[1]))

# clf = doubletdetection.BoostClassifier(n_iters=args.n_iterations, use_phenograph=pheno, standard_scaling=standard_scaling, verbose = True)
clf = doubletdetection.BoostClassifier(n_iters=args.n_iterations, standard_scaling=standard_scaling, verbose = True)
doublets = clf.fit(raw_counts).predict(p_thresh=args.p_thresh, voter_thresh=args.voter_thresh)

results = pd.Series(doublets, name="DoubletDetection_DropletType")
dataframe =pd.DataFrame(results)
dataframe.index = barcodes_df
dataframe.insert(0, 'barcodes_df', barcodes_df)
dataframe.DoubletDetection_DropletType = dataframe.DoubletDetection_DropletType.replace(1.0, "doublet")
dataframe.DoubletDetection_DropletType = dataframe.DoubletDetection_DropletType.replace(0.0, "singlet")

dataframe.to_csv(os.path.join(args.outdir,'DoubletDetection_results.txt'), sep = "\t", index = False)

### Figures ###
doubletdetection.plot.convergence(clf, save=os.path.join(args.outdir,'convergence_test.pdf'), show=False, p_thresh=args.p_thresh, voter_thresh=args.voter_thresh)

f3 = doubletdetection.plot.threshold(clf, save=os.path.join(args.outdir,'threshold_test.pdf'), show=False, p_step=6)
