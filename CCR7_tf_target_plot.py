# Aug 26 plotting CCR7 target
# Plotting tf targets of CCR7 using NicheNet database

import pandas as pd
import json
import scanpy as sc
import anndata
from PIL import Image
import matplotlib.pyplot as plt
import numpy as np 

h5_file_path = '/cluster/projects/schwartzgroup/fatema/data/V1_Human_Lymph_Node_spatial/filtered_feature_bc_matrix.h5'
positions_file_path = "/cluster/projects/schwartzgroup/fatema/data/V1_Human_Lymph_Node_spatial/spatial/tissue_positions_list.csv"
scalefactors_file_path = "/cluster/projects/schwartzgroup/fatema/data/V1_Human_Lymph_Node_spatial/spatial/scalefactors_json.json"
tissue_image_file_path = "/cluster/projects/schwartzgroup/fatema/data/V1_Human_Lymph_Node_spatial/spatial/tissue_hires_image.png"
output_image_path = "/cluster/projects/schwartzgroup/vg/spots_locations.png"
ligand_tf_matrix = pd.read_csv('/cluster/projects/schwartzgroup/vg/NEST_revision/Nichenet_dataset/ligand_tf_matrix.csv', index_col=0) # Load the ligand-TF matrix from Nichenet database

adata = sc.read_10x_h5(h5_file_path)

# Load spatial metadata
positions = pd.read_csv(positions_file_path, header=None)
positions.columns = ["barcode", "col1", "col2", "col3", "x", "y"]  # Rename columns
positions.set_index("barcode", inplace=True)
positions = positions.reindex(adata.obs_names)
adata.obsm['spatial'] = positions[["x", "y"]].values
adata.obsm['spatial'] = adata.obsm['spatial'][:, [1, 0]]  # Flip the x and y coordinates to match the NEST paper plotting


scalefactors = json.load(open("/cluster/projects/schwartzgroup/fatema/data/V1_Human_Lymph_Node_spatial/spatial/scalefactors_json.json"))
image = Image.open("/cluster/projects/schwartzgroup/fatema/data/V1_Human_Lymph_Node_spatial/spatial/tissue_hires_image.png")

# Define the spatial metadata, needed for i.e. spot size
adata.uns['spatial'] = {
    'default': {  # Use 'default' as the library ID
        'images': {
            'hires': np.array(image)
        },
        'scalefactors': scalefactors,
        'metadata': {
            'source_image_path': "/cluster/projects/schwartzgroup/fatema/data/V1_Human_Lymph_Node_spatial/spatial/tissue_hires_image.png"
        }
    }
}

# Expression normalization
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Identify transcription factors targeted by CCR7
ccr7_tfs = ligand_tf_matrix.loc['CCR7']  # each row is a TF
# all non-zero TF
ccr7_tfs = ccr7_tfs[ccr7_tfs > 0].index.tolist() 

##### Plot the expression of INDIVIDUAL gene targeted by CCR7 (there's like ~150)
for tf in ccr7_tfs: # for tf in ccr7_tfs[:10] if only plotting the first 10
    if tf in adata.var_names:
        sc.pl.spatial(adata, color=tf, title=f"Expression of {tf}", cmap="magma_r",show=False,img_key=None)
        plt.savefig(f"/cluster/projects/schwartzgroup/vg/NEST_revision//{tf}_spatial_plot.png", dpi=300)
        plt.close()

# ##### plot the summed expression of all target genes of CCR7
# # Check for duplicate gene names first
# duplicate_genes = adata.var_names[adata.var_names.duplicated()]
# print(f"Number of duplicate genes: {len(duplicate_genes)}")

# # Rename duplicate genes by appending a suffix
# adata.var_names = pd.Index([f"{gene}_{i}" if gene in duplicate_genes else gene
#                             for i, gene in enumerate(adata.var_names)])

# valid_tfs = [tf for tf in ccr7_tfs if tf in adata.var_names]
# adata.obs['CCR7_TF_sum'] = adata[:, valid_tfs].X.sum(axis=1) 

# # Plot
# sc.pl.spatial(adata, color='CCR7_TF_sum', title="Summed Expression of CCR7 TF", cmap="magma_r", show=False,img_key=None)

# plt.savefig("/cluster/projects/schwartzgroup/vg/CCR7_TF_sum.png",dpi=300)



