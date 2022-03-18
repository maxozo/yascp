#!/usr/bin/env python

__date__ = '2020-04-24'
__version__ = '0.0.1'

import argparse
import os
os.environ['NUMBA_CACHE_DIR']='/tmp'
os.environ['MPLCONFIGDIR']='/tmp'
import random
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib import cm
import scanpy as sc
import csv
import joblib  # for numpy matrix, joblib faster than pickle

from sklearn import metrics
from sklearn import preprocessing
from sklearn import model_selection
from sklearn.linear_model import LogisticRegression

# Import dask things that can be imported here.
from dask_jobqueue import LSFCluster
from dask.distributed import Client, performance_report, get_task_stream
# import dask
# import dask.array as da
# from dask_ml import preprocessing
# from dask_ml import model_selection
# import dask_ml.linear_model as dlm

# NOTE: The whole clustering pipeline could be replaced using dask as
# described here. One could do a grid search over (1) neighbors for clustering,
# (2) resolution for clustering, (3) regularization ... using the logreg score
# as the estimator to guide selection. See more below:
# https://ml.dask.org/examples/dask-glm.html

# Set seed for reproducibility
seed_value = 0
# 0. Set `PYTHONHASHSEED` environment variable at a fixed value
# os.environ['PYTHONHASHSEED']=str(seed_value)
# 1. Set `python` built-in pseudo-random generator at a fixed value
random.seed(seed_value)
# 2. Set `numpy` pseudo-random generator at a fixed value
np.random.seed(seed_value)


def start_dask_lsfcluster(cluster_size=5):
    """Start a dask cluster."""
    if cluster_size < 4:
        raise Exception('Too small of a cluster')
    # Settings for Sanger farm
    memory_in_gb = 20
    cluster = LSFCluster(
        queue='normal',
        walltime='00:30',
        log_directory='{}/dask_logs'.format(os.getcwd()),
        cores=4,
        memory='{} Gb'.format(memory_in_gb),
        mem=memory_in_gb*1e+9,  # should be in bytes
        lsf_units='mb',
        job_extra=[
            '-G team152',
            '-g /lt9/dask',
            '-R "select[mem>{}] rusage[mem={}]"'.format(
                int(memory_in_gb*1e+3),
                int(memory_in_gb*1e+3)
            )
        ],
        use_stdin=True
    )

    # View the job submission from Dask
    # cluster.job_script()

    # Scale cluster
    cluster.scale(cluster_size)

    # auto-scale between 10 and 100 jobs
    # cluster.adapt(
    #     minimum_jobs=int(cluster_size/4),
    #     maximum_jobs=cluster_size
    # )
    # cluster.adapt(maximum_memory="10 TB")  # use core/memory limits

    client = Client(
        cluster,
        timeout=120
    )
    client.wait_for_workers(n_workers=cluster_size)
    # print(client.scheduler_info()['services'])

    return cluster, client


def _create_colors(lr):
    n_cts = lr.classes_.shape[0]
    color_norm = colors.Normalize(vmin=-n_cts / 3, vmax=n_cts)
    ct_arr = np.arange(n_cts)
    ct_colors = cm.YlGnBu(color_norm(ct_arr))

    return ct_colors


def plot_roc(y_prob, y_test, lr):
    """Plot ROC curve. Based off of NaiveDE library."""
    ct_colors = _create_colors(lr)

    for i, cell_type in enumerate(lr.classes_):
        fpr, tpr, _ = metrics.roc_curve(y_test == cell_type, y_prob[:, i])
        plt.plot(fpr, tpr, c=ct_colors[i], lw=2)

    plt.plot([0, 1], [0, 1], color='k', ls=':')
    plt.xlabel('FPR')
    plt.ylabel('TPR')


def class_report(y_true, y_pred, y_score=None, average='macro'):
    """
    Build a text report showing the main classification metrics.

    Replaces sklearn.metrics.classification_report.

    Derived from:
    https://stackoverflow.com/questions/39685740/calculate-sklearn-roc-auc-score-for-multi-class
    """
    if y_true.shape != y_pred.shape:
        raise Exception(
            'Error! y_true {} is not the same shape as y_pred {}'.format(
                y_true.shape,
                y_pred.shape
            )
        )

    # Value counts of predictions
    labels, cnt = np.unique(y_pred, return_counts=True)
    n_classes = len(labels)
    pred_cnt = pd.Series(cnt, index=labels)

    # Compute precision, recall, F-measure and support for each class  The
    # precision is the ratio tp / (tp + fp) where tp is the number of true
    # positives and fp the number of false positives. The precision is
    # intuitively the ability of the classifier not to label as positive a
    # sample that is negative.  The recall is the ratio tp / (tp + fn) where tp
    # is the number of true positives and fn the number of false negatives. The
    # recall is intuitively the ability of the classifier to find all the
    # positive samples.  The F-beta score can be interpreted as a weighted
    # harmonic mean of the precision and recall, where an F-beta score reaches
    # its best value at 1 and worst score at 0.   The F-beta score weights
    # recall more than precision by a factor of beta. beta == 1.0 means recall
    # and precision are equally important.
    metrics_summary = metrics.precision_recall_fscore_support(
        y_true=y_true,
        y_pred=y_pred,
        labels=labels
    )
    avg = list(metrics.precision_recall_fscore_support(
        y_true=y_true,
        y_pred=y_pred,
        average='weighted'
    ))

    metrics_sum_index = ['precision', 'recall', 'f1-score', 'support']
    class_report_df = pd.DataFrame(
        list(metrics_summary),
        index=metrics_sum_index,
        columns=labels
    )

    support = class_report_df.loc['support']
    total = support.sum()
    class_report_df['avg / total'] = avg[:-1] + [total]

    class_report_df = class_report_df.T
    class_report_df['pred'] = pred_cnt
    class_report_df['pred'].iloc[-1] = total

    if not (y_score is None):
        fpr = dict()
        tpr = dict()
        roc_auc = dict()
        for label_it, label in enumerate(labels):
            fpr[label], tpr[label], _ = metrics.roc_curve(
                (y_true == label).astype(int),
                y_score[:, label_it]
            )
            roc_auc[label] = metrics.auc(fpr[label], tpr[label])

        # NOTE: Micro will fail if we have a classification that is missing
        # a prediction in the lm, for instance if it was not included in
        # the training set. If this happens, then
        # lb.transform(y_true).ravel() != y_score[:, 1].ravel()
        if average == 'micro':
            lb = preprocessing.LabelBinarizer()
            if len(y_true.shape) == 1:
                lb.fit(y_true)
            elif average == 'micro':
                raise Exception('Error with LabelBinarizer')

            if n_classes <= 2:
                fpr['avg / total'], tpr['avg / total'], _ = metrics.roc_curve(
                    lb.transform(y_true).ravel(),
                    y_score[:, 1].ravel()
                )
            else:
                fpr['avg / total'], tpr['avg / total'], _ = metrics.roc_curve(
                    lb.transform(y_true).ravel(),
                    y_score.ravel()
                )
            roc_auc['avg / total'] = metrics.auc(
                fpr['avg / total'],
                tpr['avg / total']
            )
        elif average == 'macro':
            # First aggregate all false positive rates
            all_fpr = np.unique(np.concatenate([
                fpr[i] for i in labels]
            ))
            # Then interpolate all ROC curves at this points
            mean_tpr = np.zeros_like(all_fpr)
            for i in labels:
                mean_tpr += np.interp(all_fpr, fpr[i], tpr[i])
            # Finally average it and compute AUC
            mean_tpr /= n_classes

            fpr['macro'] = all_fpr
            tpr['macro'] = mean_tpr

            roc_auc['avg / total'] = metrics.auc(fpr['macro'], tpr['macro'])

        class_report_df['AUC'] = pd.Series(roc_auc)

    # Catch the case where true label not predicted in lr, perhaps because
    # too few training cases.
    for i in np.unique(y_true):
        if i not in class_report_df.index:
            print(
                'Adding category ({}) from {}.'.format(
                    i,
                    'truth with no prediction to report'
                )
            )
            class_report_df = class_report_df.append(pd.Series(
                [np.nan]*len(class_report_df.columns),
                index=class_report_df.columns,
                name=i
            ))

    class_report_df = class_report_df.sort_index()

    return class_report_df


def logistic_model(
    X,
    cell_labels,
    sparsity=1.0,
    train_size_fraction=0.75,
    n_jobs=None,
    standarize=True,
    with_mean=False,
    verbose=True,
    use_dask=False,
    out_file='out_file'
):
    """Fit logistic regression model. Based off of NaiveDE."""
    # Standarize features - this is especially important for SAGA
    # https://scikit-learn.org/stable/modules/generated/sklearn.preprocessing.StandardScaler.html
    #
    # Many elements used in the objective function of a learning algorithm
    # (such as the RBF kernel of Support Vector Machines or the L1 and L2
    # regularizers of linear models) assume that all features are centered
    # around 0 and have variance in the same order. If a feature has a
    # variance that is orders of magnitude larger that others, it might
    # dominate the objective function and make the estimator unable to
    # learn from other features correctly as expected.
    #
    # NOTE: This scaler can also be applied to sparse CSR or CSC matrices by
    # passing with_mean=False to avoid breaking the sparsity structure
    # of the data.
    #
    # Standarization is also recommended by authors of
    # The Elements of Statistical Learning: Data Mining, Inference, and
    # Prediction. 9780387216065
    if standarize:
        scaler = preprocessing.StandardScaler(
            with_mean=with_mean,
            with_std=True
        )
        if with_mean and sp.sparse.issparse(X):
            X = X.todense()
        X_std = scaler.fit_transform(X)
        # X_std = (X / X.std()).dropna(1)  # Does not handle sparse matrix
        if verbose:
            print('Scaled X, with_mean={}'.format(with_mean))
    else:
        X_std = X

    # if use_dask:
    #     # Your chunks input will be normalized and stored in the third and
    #     # most explicit form. Note: chunks stands for "chunk shape" rather
    #     # than "number of chunks", so specifying chunks=1 means that you
    #     # will have many chunks, each with exactly one element.
    #     #
    #     # Automatic chunking
    #     # -1: no chunking along this dimension
    #     # None: no change to the chunking along this dimension (useful for
    #     #       rechunk)
    #     # "auto": allow the chunking in this dimension to accommodate ideal
    #     #         chunk sizes
    #     chunk_x = X.shape[0] // n_jobs
    #     chunk_y = X.shape[1] // n_jobs
    #     if verbose:
    #         print('Chunks: {},{}'.format(chunk_x, chunk_y))
    #     # NOTE: currently only allowed to chunk on x axis.
    #     X_std = da.from_array(X_std, chunks=({0: 'auto', 1: -1}))

    # Split the data into training and test data.
    # NOTE: If one does not use stratify, it is possible that a specific cell
    #       type label will be completely missing in the final dataframe.
    #       The stratify parameter makes a split so that the proportion of
    #       values in the sample produced will be the same as the proportion
    #       of values provided to parameter stratify.
    X_train, X_test, y_train, y_test = model_selection.train_test_split(
        X_std,
        cell_labels,
        stratify=cell_labels,
        random_state=61,
        train_size=train_size_fraction
    )
    # if use_dask:
    #     # Calling dask.persist will preserve our data in memory, so no
    #     # computation will be needed as we pass over our data many times.
    #     # For example if our data came from CSV files and was not persisted,
    #     # then the CSV files would have to be re-read on each pass. This is
    #     # desirable if the data does not fit in RAM, but not slows down
    #     # our computation otherwise.
    #     X_train, X_test, y_train, y_test = dask.persist(
    #         X_train, X_test, y_train, y_test
    #     )
    if verbose:
        print(
            'Split X into training {} and test {} sets.'.format(
                X_train.shape,
                X_test.shape
            )
        )

    # As noted in scanpy documentation:
    # https://icb-scanpy.readthedocs-hosted.com/en/stable/api/scanpy.tl.rank_genes_groups.html
    # penalty='l1' to try to come up with a minimal set of genes that are good
    # predictors (sparse solution meaning few non-zero fitted coefficients).
    #
    # L1 = Lasso Regression, L2 = Ridge Regression. Because Lasso shrinks the
    # less important feature’s coefficient to zero (removing some feature
    # altogether), it works well for feature selection when there are a
    # large number of features.
    # https://scikit-learn.org/stable/auto_examples/linear_model/plot_logistic_l1_l2_sparsity.html
    #
    # SAGA is a good solver for large datasets  both the number of samples
    # and the number of features) and supports l1 penalty.
    # https://stackoverflow.com/questions/38640109/logistic-regression-python-solvers-defintions
    if use_dask:
        # lr = dlm.LogisticRegression(
        #     penalty='l1',
        #     C=sparsity
        # )
        # NOTE: multi_class='ovr' may scale better than multinomial since
        # each model is fit independantly?
        lr = LogisticRegression(
            penalty='l1',
            solver='saga',  # ‘liblinear’ and ‘saga’ also handle L1 penalty
            C=sparsity,
            random_state=23,
            n_jobs=-1,
            multi_class='multinomial'
        )
        # NOTE: Could add scater = [X_train] to joblib call to give a
        # local copy of X_train to each node.
        # NOTE: If no nodes are yet available in Client, code below will throw
        # error.
        with joblib.parallel_backend('dask'):
            with performance_report(
                filename='{}-dask-ml-performance_report.html'.format(
                    out_file
                )
            ):
                with get_task_stream(
                    filename='{}-dask-ml-task_stream.html'.format(
                        out_file
                    )
                ):
                    lr.fit(X_train, y_train)
    else:
        lr = LogisticRegression(
            penalty='l1',
            solver='saga',  # ‘liblinear’ and ‘saga’ also handle L1 penalty
            C=sparsity,
            n_jobs=n_jobs,
            random_state=23,
            # max_iter=4,  # Useful for debugging.
            multi_class='multinomial'
        )
        lr.fit(X_train, y_train)
    if verbose:
        print('Completed: LogisticRegression.')

    # NOTE: P-value estimation
    # https://scikit-learn.org/stable/modules/linear_model.html#logistic-regression
    # It is possible to obtain the p-values and confidence intervals for
    # coefficients in cases of regression without penalization. The
    # statsmodels package <https://pypi.org/project/statsmodels/> natively
    # supports this. Within sklearn, one could use bootstrapping instead
    # as well.

    y_prob = lr.predict_proba(X_test)
    # if use_dask:
    #     y_prob = lr.predict_proba(X_test).compute()
    y_prob_df = pd.DataFrame(
        y_prob,
        columns=lr.classes_
    )
    y_prob_df['cell_label_true'] = y_test
    if verbose:
        print('Completed: predict_proba.')

    # Make a classifier report
    model_report = class_report(
        y_true=y_test,
        y_pred=lr.predict(X_test),
        y_score=y_prob
    )
    # Add the number of cells in each class (index) in the
    # (a) full dataset and (b) training dataset.
    categories, counts = np.unique(cell_labels, return_counts=True)
    cat_counts = dict(zip(categories, counts))
    model_report['n_cells_full_dataset'] = model_report.index.map(cat_counts)
    categories, counts = np.unique(y_train, return_counts=True)
    cat_counts = dict(zip(categories, counts))
    model_report['n_cells_training_dataset'] = model_report.index.map(
        cat_counts
    )

    # Compute p-values
    # https://github.com/pachterlab/NYMP_2018/blob/master/10x_example-logR/10x_example_logR-TCC_notebook.ipynb
    # k = 1
    # gene_score = log_loss(logr_labels,pred)
    # llf = -gene_score*(N1+N2)
    # llr = llf-llnull
    # llr_pval = stats.chi2.sf(2*llr, k) #survival function defined as 1-cdf

    # Check out the performance of the model
    score_train = lr.score(X_train, y_train)
    score_test = lr.score(X_test, y_test)
    # if use_dask:
    #     score_train = score_train.compute()
    #     score_test = score_test.compute()
    if verbose:
        print('Training score:\t{}'.format(str(score_train)))
        print('Test score:\t{}'.format(str(score_test)))

    return lr, model_report, y_prob_df, score_train, score_test


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Fits logistic regression to predict labels in adata.obs["cluster"]'
            """
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )

    parser.add_argument(
        '-h5', '--h5_anndata',
        action='store',
        dest='h5',
        required=True,
        help='H5 AnnData file where clusters have been saved to cluster slot.'
    )

    parser.add_argument(
        '-ncells', '--number_cells',
        action='store',
        dest='ncells',
        default=-1,
        type=int,
        help='Downsample to this number of cells.\
            (default: no downsampling)'
    )

    # NOTE: could potentially use a grid search estimator to set this
    # https://scikit-learn.org/stable/modules/generated/sklearn.model_selection.GridSearchCV.html#sklearn.model_selection.GridSearchCV
    parser.add_argument(
        '-s', '--sparsity',
        action='store',
        dest='sparsity',
        default=0.1,
        type=float,
        help='LogisticRegression sparsity or inverse of regularization\
            strength; must be a positive float. Like in support vector\
            machines, smaller values specify stronger regularization.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-tsf', '--train_size_fraction',
        action='store',
        dest='train_size_fraction',
        default=0.67,
        type=float,
        help='Fraction of the data to use for training set.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-tsc', '--train_size_cells',
        action='store',
        dest='train_size_cells',
        default=0,
        type=int,
        help='Number of cells to use for training set. If > 0 all\
            remaining cells not randomly selected for training will be used\
            for the test set. Overrides <train_size_fraction>.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-ncpu', '--number_cpu',
        action='store',
        dest='ncpu',
        default=1,
        type=int,
        help='Number of CPUs to use.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-d', '--dask_scale',
        action='store',
        dest='dask_scale',
        default=0,
        type=int,
        help='Scale using Dask if > 0.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='',
        help='Basename of output files, assuming output in current working \
            directory.\
            (default: <h5_anndata_file_name>)'
    )

    options = parser.parse_args()

    # Fixed settings.
    verbose = True

    # Scanpy settings
    # sc.settings.figdir = os.getcwd()  # figure output directory to match base
    # sc.settings.n_jobs = options.ncpu  # number CPUs
    # sc.settings.max_memory = 500  # in Gb
    # sc.set_figure_params(dpi_save = 300)

    # Get the out file base.
    out_file_base = options.of
    if out_file_base == '':
        out_file_base = '{}'.format(
            os.path.basename(options.h5.rstrip('h5ad')).rstrip('.')
        )

    # Load the AnnData file.
    # This file should already have clusters identified and saved to the
    # clusters slot.
    adata = sc.read_h5ad(filename=options.h5)

    # Optionally downsample cells
    if options.ncells > -1:
        n_cells_start = adata.n_obs
        if n_cells_start < options.ncells:
            raise Exception('Fewer cells than ncells specified.')
        # Downsample cells by fraction
        # sc.pp.subsample(
        #     adata,
        #     fraction=0.2,
        #     copy=False
        # )
        # Downasmple cells to a specific number.
        sc.pp.subsample(
            adata,
            n_obs=options.ncells,
            random_state=231,
            copy=False
        )
        print('Cell downsample applied: {} dropped, {} remain.'.format(
            n_cells_start - adata.n_obs,
            adata.n_obs
        ))
    categories, counts = np.unique(
        adata.obs['cluster'].values,
        return_counts=True
    )
    assert min(counts > 1), 'ERROR: smallest cluster has 1 cell.'

    # Set X to cp10k
    adata.X = np.expm1(adata.layers['log1p_cp10k'])
    # Set X to raw counts
    # adata.X = adata.layers['counts']

    # Ensure un-informative genes are filtered
    # sc.pp.filter_genes(adata, min_cells=5)

    # Subset to only highly variable genes
    subset_to_hcg = False
    if subset_to_hcg:
        adata = adata[:, adata.var['highly_variable']]

    # NOTE: We could use scanpy to stanardize the data, but just like sklearn
    #       this will also change the matrix from Compressed Sparse Row format
    #       to dense format.
    # sc.pp.scale(
    #     adata,
    #     zero_center=True,
    #     max_value=None,
    #     copy=False
    # )

    # NOTE: no need to make dense with sklearn... can keep spase
    # X = pd.DataFrame(X.toarray())

    # If train_size_cells, override the fraction so that the total number of
    # cells in the training set will be equal to train_size_cells.
    train_size_fraction = options.train_size_fraction
    if options.train_size_cells > 0:
        if options.train_size_cells >= adata.n_obs:
            raise Exception('Invalid train_size_cells.')
        train_size_fraction = (
            1 - ((adata.n_obs-options.train_size_cells)/adata.n_obs)
        )
        if verbose:
            print('Set train_size_fraction to: {}.'.format(
                train_size_fraction
            ))
    if verbose:
        print('Number cells training ({}) and testing ({}).'.format(
            int(train_size_fraction*adata.n_obs),
            int((1-train_size_fraction)*adata.n_obs)
        ))

    # Set up dask cluster if needed.
    use_dask = False
    n_jobs = options.ncpu
    if options.dask_scale > 0:
        cluster, client = start_dask_lsfcluster(options.dask_scale)
        use_dask = True
        n_jobs = -1
        print(client)

    # Here sparsity or C is the C param from
    # sklearn.linear_model.LogisticRegression
    # Inverse of regularization strength; must be a positive float.
    # Like in support vector machines, smaller values specify
    # stronger regularization.
    lr, model_report, test_results, score_train, score_test = logistic_model(
        X=adata.X,
        cell_labels=adata.obs['cluster'].values,
        sparsity=options.sparsity,
        train_size_fraction=train_size_fraction,
        n_jobs=n_jobs,
        standarize=True,
        with_mean=True,
        verbose=verbose,
        use_dask=use_dask,
        out_file=out_file_base
    )
    if use_dask:
        cluster.close()

    # Save the model
    out_f = '{}-lr_model.joblib.gz'.format(out_file_base)
    joblib.dump(
        lr,
        out_f,
        compress=('gzip', 3)
    )
    # Example of how to load a model.
    # lr = joblib.load(
    #     'CD5677E01F0E30D5-adata-normalized_pca-clustered-lr.joblib.gz'
    # )
    if verbose:
        print('Completed: save {}.'.format(out_f))

    # Save the model_report - each row cluster.
    # Add details on this call to the model report.
    model_report['sparsity'] = options.sparsity
    # Add resolution used to generate these clusters.
    for key, value in adata.uns['cluster']['params'].items():
        col_add = 'cluster__{}'.format(key)
        model_report[col_add] = value
    # Add neighbors parameters.
    for key, value in adata.uns['neighbors']['params'].items():
        col_add = 'neighbors__{}'.format(key)
        model_report[col_add] = value
    # Save.
    out_f = '{}-model_report.tsv.gz'.format(out_file_base)
    model_report.to_csv(
        out_f,
        sep='\t',
        index=True,
        index_label='class',
        quoting=csv.QUOTE_NONNUMERIC,
        na_rep='',
        compression='gzip'
    )
    if verbose:
        print('Completed: save {}.'.format(out_f))

    # Save the test results - each row is a cell and the columns are the prob
    # of that cell belonging to a particular class.
    out_f = '{}-test_result.tsv.gz'.format(out_file_base)
    test_results.to_csv(
        out_f,
        sep='\t',
        # index=True,   # NOTE: Not adding the label to test_result index.
        # index_label='cell_label',
        quoting=csv.QUOTE_NONNUMERIC,
        na_rep='',
        compression='gzip'
    )
    if verbose:
        print('Completed: save {}.'.format(out_f))

    # Plot the number of features with non-zero coefficients in each cluster.
    out_f = '{}-n_features.pdf'.format(out_file_base)
    df_plt = pd.DataFrame({
        'classes': lr.classes_,
        'features': (lr.coef_ > 0).sum(1)
    })
    df_plt = df_plt.set_index('classes')
    # Add in catgories with no predictive model (e.g., becuase they were too
    # few in training).
    for i in adata.obs['cluster'].cat.categories:
        if i not in df_plt.index:
            df_plt = df_plt.append(pd.Series(
                [0],
                index=df_plt.columns,
                name=i
            ))
    fig = plt.figure(figsize=(max(0.5*len(df_plt.index), 5), 4))
    # plt.bar(lr.classes_, n_features)
    plt.bar(df_plt.index, df_plt['features'])
    plt.xlabel('Cluster')
    plt.ylabel('Features with coefficient > 0')
    plt.xticks(rotation=90)
    for i in df_plt.index:
        plt.annotate(
            str(df_plt.loc[i, 'features']),
            xy=(i, df_plt.loc[i, 'features'])
        )
    fig.savefig(
        out_f,
        #dpi=300,
        bbox_inches='tight'
    )
    plt.close(fig)

    # Plot ROC of the test and truth.
    out_f = '{}-roc.pdf'.format(out_file_base)
    fig = plt.figure()
    cell_label_true = test_results.pop('cell_label_true')
    plot_roc(test_results.values, cell_label_true.values, lr)
    fig.savefig(
        out_f,
        #dpi=300,
        bbox_inches='tight'
    )
    plt.close(fig)
    if verbose:
        print('Completed: save {}.'.format(out_f))

    # Plot the AUC vs cluster size to see if smaller clusters have poorer AUC.
    out_f = '{}-cluster_size_auc.pdf'.format(out_file_base)
    df_plt = model_report.drop('avg / total').fillna(0)
    fig = plt.figure()
    plt.scatter(df_plt['n_cells_full_dataset'], df_plt['AUC'], alpha=0.5)
    plt.xlabel('Number of cells in cluster (full dataset)')
    plt.ylabel('AUC in training data')
    plt.ylim(0, 1)
    fig.savefig(
        out_f,
        #dpi=300,
        bbox_inches='tight'
    )
    plt.xscale('log', basex=10)
    # Add annotation of the cluster
    for index, row in df_plt.iterrows():
        plt.annotate(
            index,  # this is the text
            (row['n_cells_full_dataset'], row['AUC']),  # the point to label
            textcoords='offset points',  # how to position the text
            xytext=(0, 10),  # distance from text to points (x,y)
            ha='center'   # horizontal alignment can be left, right or center
        )
    fig.savefig(
        '{}-cluster_size_auc_log10.pdf'.format(out_file_base),
        #dpi=300,
        bbox_inches='tight'
    )
    plt.close(fig)
    if verbose:
        print('Completed: save {}.'.format(out_f))

    # Make a dataframe of the coefficients.
    # Rows = cell type label and columns = genes tested
    lr_res = pd.DataFrame.from_records(
        lr.coef_,
        index=lr.classes_,
        columns=adata.var.index
    )
    # Save the lr_res dataframe.
    out_f = '{}-lr_coef.tsv.gz'.format(out_file_base)
    lr_res.to_csv(
        out_f,
        sep='\t',
        index=True,
        index_label='cell_label',
        quoting=csv.QUOTE_NONNUMERIC,
        na_rep='',
        compression='gzip'
    )
    if verbose:
        print('Completed: save {}.'.format(out_f))


def dev():
    adata = sc.read_h5ad(filename='adata-normalized_pca-clustered.h5ad')
    # sc.pp.subsample(
    #     adata,
    #     n_obs=2000,
    #     copy=False
    # )
    categories, counts = np.unique(
        adata.obs['cluster'].values,
        return_counts=True
    )
    assert min(counts > 1), 'ERROR: smallest cluster has 1 cell.'

    adata.X = adata.layers['log1p_cp10k']

    X=adata.X
    cell_labels=adata.obs['cluster'].values
    sparsity=0.1
    train_size_fraction=0.05
    n_jobs=1
    standarize=True
    with_mean=True
    verbose=True
    use_dask=False
    out_file='test'

    lr = LogisticRegression(
        penalty='l1',
        solver='saga',  # ‘liblinear’ and ‘saga’ also handle L1 penalty
        C=sparsity,
        n_jobs=n_jobs,
        max_iter=4,
        multi_class='ovr'
    )
    lr.fit(X_train, y_train)

    y_true=y_test
    y_pred=lr.predict(X_test)
    y_score=lr.predict_proba(X_test)

    a, b = np.unique(y_pred, return_counts=True)
    print(len(b))


if __name__ == '__main__':
    main()
