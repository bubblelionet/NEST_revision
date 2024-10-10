import pandas as pd
import json
import scanpy as sc
import anndata
from PIL import Image
import matplotlib.pyplot as plt
import numpy as np
from sklearn.preprocessing import QuantileTransformer  # Import QuantileTransformer

h5_file_path = '/cluster/projects/schwartzgroup/fatema/data/V1_Human_Lymph_Node_spatial/filtered_feature_bc_matrix.h5'
positions_file_path = "/cluster/projects/schwartzgroup/fatema/data/V1_Human_Lymph_Node_spatial/spatial/tissue_positions_list.csv"
scalefactors_file_path = "/cluster/projects/schwartzgroup/fatema/data/V1_Human_Lymph_Node_spatial/spatial/scalefactors_json.json"
tissue_image_file_path = "/cluster/projects/schwartzgroup/fatema/data/V1_Human_Lymph_Node_spatial/spatial/tissue_hires_image.png"
output_image_path = "/cluster/projects/schwartzgroup/vg/spots_locations.png"

gene_to_plot = 'CCL19'

# Load AnnData object
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

# Normalization
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# # Quantile normalization on the expression matrix
# transformer = QuantileTransformer(output_distribution='uniform', random_state=42)  # or 'normal'
# adata.X = transformer.fit_transform(adata.X.toarray())  # Apply quantile normalization to expression values


# plot using scanpy
sc.pl.spatial(adata, color=gene_to_plot, size=1.5, cmap='magma_r',img_key=None) #Plot just the spots without the background image
plt.savefig("/cluster/projects/schwartzgroup/vg/spots_locations.png")



