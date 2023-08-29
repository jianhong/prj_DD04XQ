import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

import seaborn as sns

import os, sys, shutil, importlib, glob
from tqdm.notebook import tqdm

from celloracle import motif_analysis as ma
import celloracle as co

cicero_output_folder = 'cicero_output'
cicero_peak_filename = 'all_peaks.csv'
cicero_connections_filename = 'cicero_connections.csv'
ref_genome = "hg38"
integrated_coaccess = 0.8
output_filename = "processed_peak_file.csv"

peaks = pd.read_csv(os.path.join(cicero_output_folder, cicero_peak_filename),
        index_col=0)
peaks = peaks.x.values
peaks

# Load Cicero coaccessibility scores.
cicero_connections = pd.read_csv(os.path.join(cicero_output_folder, cicero_connections_filename),
        index_col=0)
cicero_connections.head()

# Load annotation
tss_annotated = ma.get_tss_info(peak_str_list=peaks, ref_genome=ref_genome)
tss_annotated.tail()

# Integrate TSS info and cicero connections
integrated = ma.integrate_tss_peak_with_cicero(tss_peak=tss_annotated,
                                               cicero_connections=cicero_connections)
print(integrated.shape)
integrated.head()

# Filter peaks
peak = integrated[integrated.coaccess >= integrated_coaccess]
peak = peak[["peak_id", "gene_short_name"]].reset_index(drop=True)

print(peak.shape)
peak.head()

peak.to_csv(output_filename)

