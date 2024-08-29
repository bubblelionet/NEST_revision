# Aug 26 plotting CCR7 target
# Plotting tf targets of CCR7 using NicheNet database

import pandas as pd
import json
import scanpy as sc
import anndata
from PIL import Image
import matplotlib.pyplot as plt
import numpy as np 

# Load the ligand-TF matrix from Nichenet database
ligand_tf_matrix = pd.read_csv('/cluster/projects/schwartzgroup/vg/NEST_revision/Nichenet_dataset/ligand_tf_matrix.csv', index_col=0)

# Identify transcription factors targeted by CCR7
ccr7_tfs = ligand_tf_matrix.loc['CCR7']
# if want to see all non-zero TF
ccr7_tfs = ccr7_tfs[ccr7_tfs > 0].index.tolist() 
# if want to see top 5 TF
ccr7_tfs = ccr7_tfs.sort_values(ascending=False)
ccr7_tfs_top15 = ccr7_tfs.head(16)
ccr7_tfs_top15 = ccr7_tfs_top15.index.tolist()


##### Plot the expression of INDIVIDUAL TF targeted by CCR7 (there's like ~150)
for tf in ccr7_tfs: # for tf in ccr7_tfs[:10] if only plotting the first 10
    if tf in adata.var_names:
        sc.pl.spatial(adata, color=tf, title=f"Expression of {tf} (CCR7 target)", cmap="viridis",show=False)
        plt.savefig(f"/cluster/projects/schwartzgroup/vg/NEST_revision//{tf}_spatial_plot.png", dpi=300,img_key=None)
        plt.close()  # Close the plot to free memory

for tf in ccr7_tfs_top15:
    if tf in adata.var_names:
        sc.pl.spatial(adata, color=tf, title=f"Expression of {tf} (CCR7 target)", cmap="viridis", show=False,img_key=None)
        plt.savefig(f"/cluster/projects/schwartzgroup/vg/NEST_revision//{tf}_spatial_plot.png", dpi=300)
        plt.close()  # Close the plot to free memory


##### To plot the sum expression of all target genes of CCR7
# Check for duplicate gene names first
duplicate_genes = adata.var_names[adata.var_names.duplicated()]
print(f"Number of duplicate genes: {len(duplicate_genes)}")

# Option 1: Remove duplicate genes (keep only the first occurrence)
adata = adata[:, ~adata.var_names.duplicated()]
# Option 2: Rename duplicate genes by appending a suffix
adata.var_names = pd.Index([f"{gene}_{i}" if gene in duplicate_genes else gene
                            for i, gene in enumerate(adata.var_names)])


# Ensure that all the TFs are present in the dataset
valid_tfs = [tf for tf in ccr7_tfs if tf in adata.var_names]

# Calculate the sum of expression for all valid CCR7 target TFs
adata.obs['CCR7_TF_sum'] = adata[:, valid_tfs].X.sum(axis=1)

# Plot the summed expression of CCR7 target TFs
sc.pl.spatial(adata, color='CCR7_TF_sum', title="Summed Expression of CCR7 Target TFs", cmap="viridis")
plt.savefig("/cluster/projects/schwartzgroup/vg/spots_locations.png",dpi=300)



