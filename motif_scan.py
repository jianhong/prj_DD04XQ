# Compare with default motifs in gimmemotifs
from gimmemotifs.motif import default_motifs
from gimmemotifs.motif import MotifConfig
from gimmemotifs.motif import read_motifs

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import seaborn as sns

import os, sys, shutil, importlib, glob
from tqdm.notebook import tqdm

import celloracle as co
from celloracle import motif_analysis as ma
from celloracle.utility import save_as_pickled_object

import argparse
parser = argparse.ArgumentParser(description="""
Author: Jianhong Ou @ duke, Aug, 2023
This source code is licensed under the MIT license
This code depends celloracle
This code will prepare the base_GRN_dataframe.parquet file""")

parser.add_argument('--processed_peaks', type=str,
                required=True,
                help='the processed peak file in csv format with columns ",peak_id,gene_short_name".')
parser.add_argument('--genome', type=str, required=True, help="The reference genome, eg hg38")
parser.add_argument('--output', type=str, default="base_GRN_dataframe.parquet", help="The output filename for final results")
parser.add_argument('--genome_dir', type=str, default="/admin/apps/community/regeneromics/conda_envs/celloracle_env/share/genomes", help="The genomes folder")
parser.add_argument('--motifs', type=str, default="JASPAR2022_vertebrates.pfm", help="The motif collection.")
parser.add_argument('--info_h5', type=str, default='test1.celloracle.tfinfo', help='The h5 output filename for tf info.')
parser.add_argument('--cores', type=int, default=1, help='number of CPUs for parallel calculation.')
parser.add_argument('--fpr', type=float, default=0.02, help='False positive rate for motif identification.')

args = parser.parse_args()


ref_genome = args.genome
processed_peak_file = args.processed_peaks
genome_dir = args.genome_dir
motif_collection = args.motifs
tfinfo_h5_output = args.info_h5
tf_parquet_output = args.output
N_CPU = args.cores
fpr = args.fpr

os.makedirs(genome_dir, exist_ok = True)

motifs =  default_motifs()

# Get folder path that stores motif data.
config = MotifConfig()
motif_dir = config.get_motif_dir()

# Get motif data names
motifs_data_name = [i for i in os.listdir(motif_dir) if i.endswith(".pfm")]
motifs_data_name.sort()
motifs_data_name

path = os.path.join(motif_dir, motif_collection)
motifs = read_motifs(path)

genome_installation = ma.is_genome_installed(ref_genome=ref_genome,
                                             genomes_dir=genome_dir)
if not genome_installation:
    import genomepy
    genomepy.install_genome(name=ref_genome, provider="UCSC", genomes_dir=genome_dir)
else:
    print(ref_genome, "is installed.")

# Load annotated peak data.
peaks = pd.read_csv(processed_peak_file, index_col=0)
peaks.head()

peaks = ma.check_peak_format(peaks, ref_genome, genomes_dir=genome_dir)

# Instantiate TFinfo object
tfi = ma.TFinfo(peak_data_frame=peaks,
                ref_genome=ref_genome,
                genomes_dir=genome_dir)

tfi.scan(fpr=fpr,
         motifs=motifs,  # If you enter None, default motifs will be loaded.
         verbose=True,
         n_cpus=N_CPU)
# Save tfinfo object
tfi.to_hdf5(file_path=tfinfo_h5_output)
# tfi=ma.load_TFinfo(tfinfo_h5_output)

# Check motif scan results
tfi.scanned_df.head()

# Reset filtering
tfi.reset_filtering()

# Do filtering
tfi.filter_motifs_by_score(threshold=10)

# Format post-filtering results.
tfi.make_TFinfo_dataframe_and_dictionary(verbose=True)

df = tfi.to_dataframe()
df.head()
# df.loc[:, 'Cebpa':'VDR'].sum(0)

# Save result as a dataframe
df.to_parquet(tf_parquet_output)
