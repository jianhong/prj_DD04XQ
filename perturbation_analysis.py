#!/usr/bin/env python

#conda activate /admin/apps/community/regeneromics/conda_envs/celloracle_env
# sample code: python ../perturbation_analysis.py --sample_id integrated --output_directory cell_oracle --min_mass 2.4 --scale 20
# python ../perturbation_analysis.py --sample_id integrated --output_directory cell_oracle_BRPC21091 --min_mass 2.4 --scale 20
import os
import sys
import argparse

import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

import celloracle as co

parser = argparse.ArgumentParser(description="""
Author: Jianhong Ou @ duke, Aug, 2023
This source code is licensed under the MIT license
This code depends celloracle
This code will do GRN model construction and network analysis""")

parser.add_argument('--sample_id', type=str, required=True, help="The cell ranger output folder")
parser.add_argument('--output_directory', type=str, default="cell_oracle", help="The output folder")
parser.add_argument('--min_mass', type=float, default=0, help="The min mass for plot")
parser.add_argument('--scale', type=float, default=0, help="The scale for the arrow length")

args = parser.parse_args()

output_directory = args.output_directory
sample_id = args.sample_id
scale_simulation = args.scale
min_mass = args.min_mass
save_folder = os.path.join(output_directory, "perturbation_analysis_figures")
os.makedirs(save_folder, exist_ok=True)
oracle = co.load_hdf5(os.path.join(output_directory, sample_id+'.celloracle.oracle'))
links = co.load_hdf5(os.path.join(output_directory, "links.celloracle.links"))

links.filter_links()
oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
oracle.fit_GRN_for_simulation(alpha=10,
                              use_cluster_specific_TFdict=True)

# Check gene expression
sorted_links_score_top = links.merged_score[links.merged_score['eigenvector_centrality']==links.merged_score['eigenvector_centrality'].max()]
top_factors = list(set(list(sorted_links_score_top.index)))
os.makedirs(os.path.join(save_folder, "scale"), exist_ok = True)
for goi in top_factors:
    # get 90% expression level for overExpression analysis
    per90 = np.percentile(sc.get.obs_df(oracle.adata, keys=[goi], layer="imputed_count").values, 90)
    for expr,epr_level in {'KO':0.0, 'OE':per90}.items():
        sc.pl.draw_graph(oracle.adata, color=[goi, oracle.cluster_column_name],
                         layer="imputed_count", use_raw=False, cmap="viridis")
        sc.get.obs_df(oracle.adata, keys=[goi], layer="imputed_count").hist()
        plt.savefig(os.path.join(save_folder, "histogram."+goi+".pdf"))
        # Enter perturbation conditions to simulate signal propagation after the perturbation.
        try:
            oracle.simulate_shift(perturb_condition={goi: epr_level},
                                  n_propagation=3)
        except:
            continue
        # check the simulation value
        oracle.evaluate_and_plot_simulation_value_distribution(
            n_genes=4,
            save=os.path.join(save_folder, f"{expr}_evaluate_{goi}_simulation_value_distribution_top_4_genes"))
        # Get transition probability
        oracle.estimate_transition_prob(n_neighbors=200,
                                        knn_random=True,
                                        sampled_fraction=1)
        # Calculate embedding
        oracle.calculate_embedding_shift(sigma_corr=0.05)
        if goi == top_factors[0]:
            for scale in [50, 25, 20, 15, 10, 5, 2, 1]:
                fig, ax = plt.subplots(1, 2,  figsize=[26, 12])
                # Show quiver plot
                oracle.plot_quiver(scale=scale, ax=ax[0])
                ax[0].set_title(f"Simulated cell identity shift vector: {goi} {expr}")
                # Show quiver plot that was calculated with randomized graph.
                oracle.plot_quiver_random(scale=scale, ax=ax[1])
                ax[1].set_title(f"Randomized simulation vector")
                plt.savefig(os.path.join(save_folder, "scale", f"{expr}_quiver_scale_"+str(scale)+".pdf"))
        # n_grid = 40 is a good starting value.
        n_grid = 40
        oracle.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
        
        if goi == top_factors[0]:
            # Search for best min_mass.
            oracle.suggest_mass_thresholds(n_suggestion=12)
            plt.savefig(os.path.join(save_folder, f"{expr}_suggest_mass_thresholds.pdf"))
            # According to the results, find the optimal min_mass which can cover all the cluster
            if min_mass == 0:
                min_mass = float(input(f'Please input the minimal mass value by the {save_folder}/{expr}_suggest_mass_thresholds.pdf: '))
            oracle.calculate_mass_filter(min_mass=min_mass, plot=True)
            if scale_simulation == 0:
                while 1:
                    new_scale = float(input(f"Please input the simulate scale for the quiver plot \nInput 0 to quit fix the scale: "))
                    if new_scale == 0.0:
                        break
                    scale_simulation = new_scale
                    fig, ax = plt.subplots(1, 2,  figsize=[26, 12])
                    # Show quiver plot
                    oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax[0])
                    ax[0].set_title(f"Simulated cell identity shift vector: {goi} {expr}")
                    # Show quiver plot that was calculated with randomized graph.
                    oracle.plot_simulation_flow_random_on_grid(scale=scale_simulation, ax=ax[1])
                    ax[1].set_title(f"Randomized simulation vector")
                    plt.savefig(os.path.join(save_folder, f"{expr}_quiver_plot_{goi}_{min_mass}_{scale_simulation}.pdf"))
            else:
                fig, ax = plt.subplots(1, 2,  figsize=[26, 12])
                # Show quiver plot
                oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax[0])
                ax[0].set_title(f"Simulated cell identity shift vector: {goi} {expr}")
                # Show quiver plot that was calculated with randomized graph.
                oracle.plot_simulation_flow_random_on_grid(scale=scale_simulation, ax=ax[1])
                ax[1].set_title(f"Randomized simulation vector")
                plt.savefig(os.path.join(save_folder, f"{expr}_quiver_plot_{goi}_{min_mass}_{scale_simulation}.pdf"))
        else:
                fig, ax = plt.subplots(1, 2,  figsize=[26, 12])
                # Show quiver plot
                oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax[0])
                ax[0].set_title(f"Simulated cell identity shift vector: {goi} {expr}")
                # Show quiver plot that was calculated with randomized graph.
                oracle.plot_simulation_flow_random_on_grid(scale=scale_simulation, ax=ax[1])
                ax[1].set_title(f"Randomized simulation vector")
                plt.savefig(os.path.join(save_folder, f"{expr}_quiver_plot_{goi}_{min_mass}_{scale_simulation}.pdf"))
        # Plot vector field with cell cluster
        fig, ax = plt.subplots(figsize=[12, 12])
        
        oracle.plot_cluster_whole(ax=ax, s=10)
        oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax, show_background=False)
        plt.savefig(os.path.join(save_folder, f"{expr}_quiver_plot_{goi}_{min_mass}_{scale_simulation}_cluster.pdf"))



