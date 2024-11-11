# Aug 25 Plotting Lymph Visium data 

import pandas as pd
import json
import scanpy as sc
import anndata
from PIL import Image
import matplotlib.pyplot as plt
import numpy as np

# Load AnnData object (if not already loaded)
adata = sc.read_10x_h5('/cluster/projects/schwartzgroup/fatema/data/V1_Human_Lymph_Node_spatial/filtered_feature_bc_matrix.h5')

# Load spatial metadata
positions = pd.read_csv("/cluster/projects/schwartzgroup/fatema/data/V1_Human_Lymph_Node_spatial/spatial/tissue_positions_list.csv", header=None)
positions.columns = ["barcode", "col1", "col2", "col3", "x", "y"]  # Rename columns as needed
positions.set_index("barcode", inplace=True)
positions = positions.reindex(adata.obs_names)
adata.obsm['spatial'] = positions[["x", "y"]].values

# Flip the x and y coordinates to match the NEST paper plotting
adata.obsm['spatial'] = adata.obsm['spatial'][:, [1, 0]]

scalefactors = json.load(open("/cluster/projects/schwartzgroup/fatema/data/V1_Human_Lymph_Node_spatial/spatial/scalefactors_json.json"))
image = Image.open("/cluster/projects/schwartzgroup/fatema/data/V1_Human_Lymph_Node_spatial/spatial/tissue_hires_image.png")

### This step is not neccesary for plotting ###
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



# plot using scanpy
sc.pl.spatial(adata, color='CCL19', size=1.5, img_key=None) #Plot just the spots without the background image
plt.savefig("/cluster/projects/schwartzgroup/vg/spots_locations.png")



