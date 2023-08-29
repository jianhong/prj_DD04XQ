#!/usr/bin/env python

#conda activate /admin/apps/community/regeneromics/conda_envs/celloracle_env

import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib import rcParams
import scanpy as sc
import os
import argparse

sc.settings.verbosity = 3  # verbosity: errors (0), warnings (1), info (2), hints (3)
# sc.logging.print_versions()

parser = argparse.ArgumentParser(description="""
Author: Jianhong Ou @ duke, Aug, 2023
This source code is licensed under the MIT license
This code depends scanpy
This code will preprocess the scRNAseq data""")

parser.add_argument('--sample_id',
                type=str, required=True,
                help="The cell ranger output folder")
parser.add_argument('--infile', type=str,
                required=True,
                help='The input h5ad file.')
parser.add_argument('--output_sub_directory',
                type=str, default="cell_oracle",
                help="The output folder")
parser.add_argument('--output_file_postfix',
                type=str, default="processed.h5ad",
                help="The output postfix of processed h5ad")
parser.add_argument('--cutoff_gene_number',
                type=int, default=7000,
                help="filter condition for gene number")
parser.add_argument('--cutoff_mt',
                type=int, default=20,
                help="filter condition for mitochondrial gene percentage")
parser.add_argument('--mt_prefix',
                type=str, default='MT-',
                help="prefix for mitochondrial gene")
parser.add_argument('--n_neighbors',
                type=int, default=10,
                help="number of neighbors")
parser.add_argument('--n_pcs',
                type=int, default=40,
                help="number of pcs")
parser.add_argument('--resolution',
                type=float, default=0.8,
                help="cluster resolution")


args = parser.parse_args()

output_directory = os.path.join(args.sample_id, args.output_sub_directory)
results_file = args.infile
output_file = args.sample_id+'.'+args.output_file_postfix
cutoff_gene = args.cutoff_gene_number
cutoff_mt = args.cutoff_mt
n_neighbors= args.n_neighbors
n_pcs= args.n_pcs
resolution= args.resolution
gene_markers = ['SFTPC','SFTPB','SCGB3A2','KRT5','FOXJ1','AGER','LTF','EPCAM','MUC5B','SCGB1A1','CLDN18']
sc.settings.figdir = output_directory
adata = sc.read_h5ad(results_file)
## Preprocessing
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
adata.var['mt'] = adata.var_names.str.startswith(args.mt_prefix)  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pl.violin(adata, [ 'pct_counts_mt'],
             jitter=0.4, multi_panel=True,
             save="_pct_mt.pdf")
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt',
            save="_total_counts.vs.pct_mt.pdf")
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts',
            save="_total_counts.vs.gene_counts.pdf")
adata = adata[adata.obs.n_genes_by_counts < cutoff_gene, :]
adata = adata[adata.obs.pct_counts_mt < cutoff_mt, :]
sc.pp.normalize_total(adata, target_sum=1e4)
adata.raw = adata
adata.layers["raw_count"] = adata.raw.X.copy()
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata, save="_highy_variable_genes.pdf")
adata = adata[:, adata.var.highly_variable]
sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca(adata, color='SFTPC', save="_pca.pdf")
sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
sc.tl.diffmap(adata)
sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep='X_diffmap')
sc.tl.umap(adata)
sc.tl.leiden(adata)
sc.pl.umap(adata, color=['leiden'],legend_loc='on data', save="_leiden.pdf")
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False,
        save='_top25_gene_by_ttest.pdf')
sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False,
        save='_top25_gene_by_wilcoxon.pdf')
sc.tl.louvain(adata, resolution=resolution)
sc.tl.paga(adata, groups='louvain')
sc.pl.paga(adata, color='louvain', save="_louvain.pdf")
sc.tl.draw_graph(adata, init_pos='paga', random_state=123)
sc.pl.draw_graph(adata, color='louvain', legend_loc='on data', save="_louvain.pdf")
sc.pl.umap(adata, color=gene_markers,
    save="_markers.pdf")
adata.write(os.path.join(output_directory, output_file))

            

