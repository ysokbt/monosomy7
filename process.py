#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 25 22:15:09 2020

@author: drmz
"""


import numpy as np
import pandas as pd
import seaborn as sns
import umap
import proplot as plot
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib as mpl
from mpl_toolkits import mplot3d

import scanpy as sc

adata = sc.read_10x_mtx('Research/mn_scrnaseq/data/')

# Adata
sc.pl.highest_expr_genes(adata, n_top=20)
sc.pp.calculate_qc_metrics(adata, percent_top=None, log1p=False, inplace=True)
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')
sc.pp.filter_genes(adata, min_cells=100)
sc.pp.filter_cells(adata, min_genes=800)
adata.var['MT'] = adata.var.index.str.startswith('MT-') 
sc.pp.calculate_qc_metrics(adata, qc_vars=['MT'], percent_top=None, log1p=False, inplace=True)
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts')
sc.pl.scatter(adata, x='total_counts', y='pct_counts_MT')

np.count_nonzero(adata.obs['pct_counts_MT'] < 10)
adata = adata[adata.obs.pct_counts_MT < 10, :]
adata = adata[adata.obs.total_counts < 50000, :]
adata = adata[adata.obs.total_counts > 5000, :]
adata.shape


# Save Filtered Counts
norm_mat = adata.X.toarray()
feat_ft = adata.var
obs_ft = adata.obs

feat_ft.to_csv('Research/mn_scrnaseq/data/expr_feat_full.tsv.gz', index=True, sep='\t')
obs_ft.to_csv('Research/mn_scrnaseq/data/obs_full.tsv.gz', index=True, sep='\t')
np.savetxt('Research/mn_scrnaseq/data/count_mat_full.tsv.gz', X=norm_mat, delimiter='\t')

sc.pp.normalize_total(adata, target_sum=None, exclude_highly_expressed=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_disp=0.1, n_bins=50)
sc.pl.highly_variable_genes(adata)

np.count_nonzero(adata.var.highly_variable)
adata = adata[:, adata.var.highly_variable]

sc.pp.neighbors(adata, n_neighbors=30, method='umap', n_pcs=None, use_rep='X')
sc.tl.leiden(adata)
sc.tl.paga(adata, groups='leiden')
sc.pl.paga(adata, color='leiden')
sc.tl.umap(adata, n_components=2, init_pos='paga')
sc.pl.umap(adata, color='leiden')

norm_mat = adata.X
feat_ft = adata.var
obs_ft = adata.obs
umap_ft = adata.obsm['X_umap']

feat_ft.to_csv('Research/mn_scrnaseq/data/expr_feat.tsv.gz', index=True, sep='\t')
obs_ft.to_csv('Research/mn_scrnaseq/data/obs_ft.tsv.gz', index=True, sep='\t')
np.savetxt('Research/mn_scrnaseq/data/expr_mat.tsv.gz', X=norm_mat.toarray(), delimiter='\t')
np.savetxt('Research/mn_scrnaseq/data/umap_ft.tsv.gz', X=umap_ft, delimiter='\t')
