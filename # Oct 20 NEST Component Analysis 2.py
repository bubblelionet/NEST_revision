# Plotting Component Specific things on Lymoph Node

import pandas as pd
import json
import scanpy as sc
import anndata
from PIL import Image
import matplotlib.pyplot as plt
import numpy as np

# Load AnnData object
adata = sc.read_10x_h5('/cluster/projects/schwartzgroup/fatema/data/V1_Human_Lymph_Node_spatial/filtered_feature_bc_matrix.h5')
positions = pd.read_csv("/cluster/projects/schwartzgroup/fatema/data/V1_Human_Lymph_Node_spatial/spatial/tissue_positions_list.csv", header=None)
positions.columns = ["barcode", "col1", "col2", "col3", "x", "y"]  # Rename columns as needed
positions.set_index("barcode", inplace=True)
positions = positions.reindex(adata.obs_names)
adata.obsm['spatial'] = positions[["x", "y"]].values
adata.obsm['spatial'] = adata.obsm['spatial'][:, [1, 0]] # Flip the x and y coordinates to match the NEST paper plotting
scalefactors = json.load(open("/cluster/projects/schwartzgroup/fatema/data/V1_Human_Lymph_Node_spatial/spatial/scalefactors_json.json"))
image = Image.open("/cluster/projects/schwartzgroup/fatema/data/V1_Human_Lymph_Node_spatial/spatial/tissue_hires_image.png")
# This step is not neccesary for plotting
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

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
# Load NEST output csv
output_csv = pd.read_csv('/cluster/projects/schwartzgroup/vg/NEST_revision/NEST_lymph_output/V1_Human_Lymph_Node_spatial_ccc_list_top20p_TCell.csv')
# Fileter output_csv by component 5
component_5_df = output_csv[output_csv['component'] == 5]
component_5_barcodes = pd.concat([component_5_df['from_cell'], component_5_df['to_cell']]).unique()
# Filter output_csv by component 9
component_9_df = output_csv[output_csv['component'] == 9]
component_9_barcodes = pd.concat([component_9_df['from_cell'], component_9_df['to_cell']]).unique()
# Filter output_csv by component 14
component_14_df = output_csv[output_csv['component'] == 14]
component_14_barcodes = pd.concat([component_14_df['from_cell'], component_14_df['to_cell']]).unique()

# Create groups in anndata
adata.obs['Component'] = 'Other'  # Default group is 'Other'
# Update group labels for barcodes present in component_5_barcodes, component_9_barcodes, and component_14_barcodes
adata.obs.loc[adata.obs.index.isin(component_5_barcodes), 'Component'] = '5'
adata.obs.loc[adata.obs.index.isin(component_9_barcodes), 'Component'] = '9'
adata.obs.loc[adata.obs.index.isin(component_14_barcodes), 'Component'] = '14'
print(adata.obs['Component'].value_counts())

# DEG analysis
# sc.tl.rank_genes_groups(adata, groupby='Component', reference='Other', method='wilcoxon')
sc.pp.log1p(adata)  # This will transform `adata.X` to log(X + 1)
sc.pp.normalize_total(adata, target_sum=1e4)


sc.tl.rank_genes_groups(adata, groupby='Component', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)

# Get the results for each comparison
deg_results = adata.uns['rank_genes_groups']
deg_df = pd.DataFrame({group: deg_results['names'][group] for group in deg_results['names'].dtype.names})
# Add the log fold changes and p-values, add columns for each group
full_deg_df = pd.concat(
    [
        deg_df,
        pd.DataFrame({group: deg_results['logfoldchanges'][group] for group in deg_results['logfoldchanges'].dtype.names}),
        pd.DataFrame({group: deg_results['pvals'][group] for group in deg_results['pvals'].dtype.names}),
    ],
    axis=1,
)
# Save the results
deg_df.to_csv('/cluster/projects/schwartzgroup/vg/NEST_revision/LUAD_component_deg_results.csv')