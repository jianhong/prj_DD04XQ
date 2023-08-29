#!/usr/bin/env python

#conda activate /admin/apps/community/regeneromics/conda_envs/celloracle_env

import pandas as pd
import numpy as np
import scanpy as sc
import warnings
import anndata
import os
import argparse

warnings.filterwarnings('ignore')
#sc.logging.print_header()

parser = argparse.ArgumentParser(description="""
Author: Jianhong Ou @ duke, Aug, 2023
This source code is licensed under the MIT license
This code depends scanpy
This code will prepare the h5ad file for scranpy analysis""")

parser.add_argument('--filtered_matrix_h5', type=str,
                default="outs/filtered_feature_bc_matrix.h5",
                help='the cell range outs of filtered feature matrix in h5 format. It should be saved in sample_id folder')
parser.add_argument('--sample_id', type=str, required=True, help="The cell ranger output folder")
parser.add_argument('--output_sub_directory', type=str, default="cell_oracle", help="The output folder")

args = parser.parse_args()

filtered_matrix_h5 = os.path.join(args.sample_id, args.filtered_matrix_h5)
output_directory = os.path.join(args.sample_id, args.output_sub_directory)
sample_id = args.sample_id
os.makedirs(output_directory)

adata = sc.read_10x_h5(filtered_matrix_h5)
adata.obs['sample_id'] = sample_id
adata.var_names_make_unique()
#adata

adata.write(os.path.join(output_directory, sample_id + '.h5ad'), compression = 'gzip')

