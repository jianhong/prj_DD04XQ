#!/usr/bin/env python

#conda activate /admin/apps/community/regeneromics/conda_envs/celloracle_env

import os
import sys

import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import argparse


import celloracle as co
# co.__version__
#'0.14.0'

parser = argparse.ArgumentParser(description="""
Author: Jianhong Ou @ duke, Aug, 2023
This source code is licensed under the MIT license
This code depends celloracle
This code will do GRN model construction and network analysis""")

parser.add_argument('--processed_h5ad', type=str,
                required=True,
                help='the cell range outs of filtered feature matrix in h5 format.')
parser.add_argument('--sample_id', type=str, required=True, help="The cell ranger output folder")
parser.add_argument('--output_directory', type=str, default="cell_oracle", help="The output folder")
parser.add_argument('--cluster_column_name', type=str, default="louvain", help="The cluster column name")
parser.add_argument('--embedding_name', type=str, default="X_draw_graph_fa", help="The embeding name")
parser.add_argument('--GRN', type=str, default='human', help='The cluster-specific GRNs. "human" or "mouse" will load default human/mouse GRN. see https://morris-lab.github.io/CellOracle.documentation/modules/celloracle.data.html for more options.')

args = parser.parse_args()

output_directory = args.output_directory
results_file = args.processed_h5ad
output_file = args.sample_id+'.celloracle.oracle'
cluster_column_name=args.cluster_column_name
embedding_name=args.embedding_name
if args.GRN == 'human':
    base_GRN = co.data.load_human_promoter_base_GRN() #human data
elif args.GRN == 'mouse':
    base_GRN = co.data.load_mouse_scATAC_atlas_base_GRN() #mouse data
elif args.GRN == 'Celegans':
    base_GRN = co.data.load_Celegans_promoter_base_GRN()
elif args.GRN == 'Pig':
    base_GRN = co.load_Pig_promoter_base_GRN()
elif args.GRN == 'Scerevisiae':
    base_GRN = co.load_Scerevisiae_promoter_base_GRN()
elif args.GRN == 'chicken':
    base_GRN = co.load_chicken_promoter_base_GRN()
elif args.GRN == 'drosophila':
    base_GRN = co.load_drosophila_promoter_base_GRN()
elif args.GRN == 'rat':
    base_GRN = co.load_rat_promoter_base_GRN()
elif args.GRN == 'xenopus_tropicalis':
    base_GRN = co.load_xenopus_tropicalis_promoter_base_GRN()
elif args.GRN == 'zebrafish':
    base_GRN = co.load_zebrafish_promoter_base_GRN()
else:
    base_GRN = pd.read_parquet(args.GRN)

os.makedirs(output_directory, exist_ok=True)
sc.settings.figdir = output_directory
adata = sc.read_h5ad(results_file)
oracle = co.Oracle()
adata.X = adata.layers["raw_count"].copy()
# Instantiate Oracle object.
oracle.import_anndata_as_raw_count(adata=adata,
                                   cluster_column_name=cluster_column_name,
                                   embedding_name=embedding_name)
oracle.import_TF_data(TF_info_matrix=base_GRN)
#oracle

# Perform PCA
oracle.perform_PCA()

# Select important PCs
plt.plot(np.cumsum(oracle.pca.explained_variance_ratio_)[:100])
n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))>0.002))[0][0]
plt.axvline(n_comps, c="k")
plt.savefig(os.path.join(output_directory, "n_comps.pdf"))
print(n_comps)
n_comps = min(n_comps, 50)
n_cell = oracle.adata.shape[0]
k = int(0.025*n_cell)
oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8,
                      b_maxl=k*4, n_jobs=4)
oracle.to_hdf5(os.path.join(output_directory, output_file))
### Load file.
#oracle = co.load_hdf5("celloracle.oracle")

sc.pl.draw_graph(oracle.adata, color=cluster_column_name,
                save="_oracle_clustering.pdf")
# Calculate GRN for each population in "louvain_annot" clustering unit.
# This step may take some time.(~30 minutes)
links = oracle.get_links(cluster_name_for_GRN_unit=cluster_column_name, alpha=10,
                         verbose_level=10)

# Save as csv
os.makedirs(os.path.join(output_directory, "raw_GRN_for_cluster"))
for cluster in links.links_dict.keys():
    links.links_dict[cluster].to_csv(os.path.join(output_directory, "raw_GRN_for_cluster", f"raw_GRN_for_{cluster}.csv"))

links.filter_links(p=0.001, weight="coef_abs", threshold_number=2000)
links.plot_degree_distributions(plot_model=True,
                                save=f"{output_directory}/degree_distribution/")
# Calculate network scores.
links.get_network_score()
links.merged_score.head()
links.merged_score.shape
links.merged_score.to_csv(os.path.join(output_directory, "links.merged_score.csv"))
links.to_hdf5(file_path=os.path.join(output_directory, "links.celloracle.links"))
sorted_links_score = links.merged_score.sort_values(by=['eigenvector_centrality'], ascending=False)
sorted_links_score_top = sorted_links_score[sorted_links_score['eigenvector_centrality']==sorted_links_score['eigenvector_centrality'].max()]
top_factors = set(list(sorted_links_score_top.index))
# top 10 in eigenvector_centrality
for factor in top_factors:
    links.plot_score_per_cluster(goi=factor, save=f"{output_directory}/network_score_per_gene/")
# Plot degree_centrality
links.plot_score_discributions(values=["degree_centrality_all"],
                               method="boxplot",
                               save=f"{output_directory}",
                              )
links.plot_score_discributions(values=["eigenvector_centrality"],
                               method="boxplot",
                               save=f"{output_directory}")
links.plot_network_entropy_distributions(save=f"{output_directory}")

# Visualize top n-th genes with high scores.
os.makedirs(os.path.join(output_directory, "link_top30"))
for cluster in links.links_dict.keys():
    links.plot_scores_as_rank(cluster=cluster, n_gene=30, save=os.path.join(output_directory, "link_top30", f"{cluster}_ranked_score"))

# Compare GRN score between two clusters
for cluster1 in links.links_dict.keys():
    for cluster2 in links.links_dict.keys():
        links.plot_score_comparison_2D(value="eigenvector_centrality",
                               cluster1=cluster1, cluster2=cluster2,
                               percentile=98,
                               save=f"{output_directory}/score_comparison")
        links.plot_score_comparison_2D(value="betweenness_centrality",
                               cluster1=cluster1, cluster2=cluster2,
                               percentile=98,
                               save=f"{output_directory}/score_comparison")
        links.plot_score_comparison_2D(value="degree_centrality_all",
                               cluster1=cluster1, cluster2=cluster2,
                               percentile=98,
                               save=f"{output_directory}/score_comparison")


