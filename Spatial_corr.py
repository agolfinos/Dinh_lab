# -*- coding: utf-8 -*-

# load packages

import pandas as pd
import sys
import scanpy as sc
import warnings
import spatialcorr

sys.path.append('../../spatialcorr')
warnings.simplefilter(action='ignore', category=FutureWarning)


# set variables

path = '/Volumes/hqdinh2/Projects/HNC_SPORE/SpatialCorr/'

ck17_5 = sc.read_h5ad(filename = path + 'ck17_5_2023-01-04.h5ad')
ck17_7 = sc.read_h5ad(filename = path + 'ck17_7_2023-01-04.h5ad')
ck17_21 = sc.read_h5ad(filename = path + 'ck17_21_2023-01-04.h5ad')
ck17_27 = sc.read_h5ad(filename = path + 'ck17_27_2023-01-04.h5ad')

# reading in the CellChat database

df = pd.read_csv(path + 'refined_cellchat_db.csv')


# Run kernel diagnostics

spatialcorr.wrappers.kernel_diagnostics(
    ck17_5,
    'BayesSpace',
    bandwidth=5,
    contrib_thresh=10,
    dsize=20,
    fpath= path + 'ck17_5_kernel_diagnostic.png',
    fformat='png',
    dpi=150
)

spatialcorr.wrappers.kernel_diagnostics(
    ck17_7,
    'BayesSpace',
    bandwidth=5,
    contrib_thresh=10,
    dsize=20,
    fpath= path + 'ck17_7_kernel_diagnostic.png',
    fformat='png',
    dpi=150
)

spatialcorr.wrappers.kernel_diagnostics(
    ck17_21,
    'BayesSpace',
    bandwidth=5,
    contrib_thresh=10,
    dsize=20,
    fpath= path + 'ck17_21_kernel_diagnostic.png',
    fformat='png',
    dpi=150
)
spatialcorr.wrappers.kernel_diagnostics(
    ck17_27,
    'BayesSpace',
    bandwidth=5,
    contrib_thresh=10,
    dsize=20,
    fpath= path + 'ck17_27_kernel_diagnostic.png',
    fformat='png',
    dpi=150
)


for index, row in df.iterrows():
    #GOI = [row['l'], row['r']]
    try: 
        fig_name = '/Volumes/hqdinh2/Projects/HNC_SPORE/SpatialCorr/' + row['l'] + '_' + row['r'] + '_analysis_pipeline_set.png'
        spatialcorr.analysis_pipeline_pair(adata = ck17_5, gene_1 = row['l'], gene_2 = row['r'], bandwidth = 5,
        cond_key ='BayesSpace', max_perms=50, dsize=20, verbose=0, fig_path = fig_name, fig_format = 'png', dpi = 300)
    except: 
        pass

for index, row in df.iterrows():
    #GOI = [row['l'], row['r']]
    try: 
        fig_name = '/Volumes/hqdinh2/Projects/HNC_SPORE/SpatialCorr/' + row['l'] + '_' + row['r'] + '_analysis_pipeline_set.png'
        spatialcorr.analysis_pipeline_pair(adata = ck17_7, gene_1 = row['l'], gene_2 = row['r'], bandwidth = 5,
        cond_key ='BayesSpace', max_perms=50, dsize=20, verbose=0, fig_path = fig_name, fig_format = 'png', dpi = 300)
    except: 
        pass

for index, row in df.iterrows():
    #GOI = [row['l'], row['r']]
    try: 
        fig_name = '/Volumes/hqdinh2/Projects/HNC_SPORE/SpatialCorr/' + row['l'] + '_' + row['r'] + '_analysis_pipeline_set.png'
        spatialcorr.analysis_pipeline_pair(adata = ck17_21, gene_1 = row['l'], gene_2 = row['r'], bandwidth = 5,
        cond_key ='BayesSpace', max_perms=50, dsize=20, verbose=0, fig_path = fig_name, fig_format = 'png', dpi = 300)
    except: 
        pass

for index, row in df.iterrows():
    #GOI = [row['l'], row['r']]
    try: 
        fig_name = '/Volumes/hqdinh2/Projects/HNC_SPORE/SpatialCorr/' + row['l'] + '_' + row['r'] + '_analysis_pipeline_set.png'
        spatialcorr.analysis_pipeline_pair(adata = ck17_27, gene_1 = row['l'], gene_2 = row['r'], bandwidth = 5,
        cond_key ='BayesSpace', max_perms=50, dsize=20, verbose=0, fig_path = fig_name, fig_format = 'png', dpi = 300)
    except: 
        pass