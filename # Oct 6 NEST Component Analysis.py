# Plotting Component Specific things on LUAD
# To show that the target of CCC predicted that belong to a specific component is indeed upregulated

import pandas as pd
import json
import scanpy as sc
import anndata
from PIL import Image
import matplotlib.pyplot as plt
import numpy as np

# Load AnnData object
LUAD_data_folder = '/cluster/projects/schwartzgroup/fatema/data/LUAD/LUAD_GSM5702473_TD1'
adata = sc.read_10x_mtx(LUAD_data_folder, var_names='gene_symbols', make_unique=True)
positions = pd.read_csv("/cluster/projects/schwartzgroup/fatema/data/LUAD/LUAD_GSM5702473_TD1/GSM5702473_TD1_tissue_positions_list.csv", header=None)
positions.columns = ["barcode", "col1", "col2", "col3", "x", "y"]  # Rename columns as needed
positions.set_index("barcode", inplace=True)
positions = positions.reindex(adata.obs_names)
adata.obsm['spatial'] = positions[["x", "y"]].values 
adata.obsm['spatial'] = adata.obsm['spatial'][:, [1, 0]] # flip x and y


# Load the component specific target genes
# ligand_tf_matrix = pd.read_csv('/cluster/projects/schwartzgroup/vg/NEST_revision/Nichenet_dataset/ligand_tf_matrix.csv', index_col=0)
# Load NEST output csv
output_csv = pd.read_csv('/cluster/projects/schwartzgroup/fatema/NEST/output/LUAD_TD1/LUAD_TD1_ccc_list_top5000.csv')

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

######### Pairwise comparisons
sc.tl.rank_genes_groups(adata, groupby='Component', reference='9', groups=['5'], method='wilcoxon')
sc.tl.rank_genes_groups(adata, groupby='Component', reference='14', groups=['5'], method='wilcoxon')
sc.tl.rank_genes_groups(adata, groupby='Component', reference='14', groups=['9'], method='wilcoxon')

# Plot top 10 genes for component 5 on tissue
top_genes_component_5 = adata.uns['rank_genes_groups']['names']['9'][:50]  # Get the top 10 genes
print("Top 10 genes for Component 5:", top_genes_component_5)

sc.pl.spatial(adata, color=top_genes_component_5, cmap='magma_r',  ncols=5, spot_size=100, size=1.5, legend_loc='on data',img_key=None)
# sc.pl.spatial(adata, color='CCR7_TF_sum', title="Summed Expression of CCR7 Downstream Targets", cmap="magma", show=False,img_key=None)
# Get the current figure and axis
# fig = plt.gcf()
# ax = plt.gca()
# # Set the background color to dark gray
# ax.set_facecolor('darkgray')
# Save the plot
plt.savefig('/cluster/projects/schwartzgroup/vg/NEST_revision/component_9_top50_genes.png', dpi=800)


####### Want to just plot one component on the tissue #########
highlight_component = '14'
adata.obs['Component_highlight'] = adata.obs['Component'].copy()
adata.obs['Component_highlight'] = adata.obs['Component_highlight'].where(
    adata.obs['Component_highlight'] == highlight_component, 'Other'
)
palette = {
    'Other': 'darkgray',
    '5': 'darkgray',   # Color for Component 5
    '9': 'darkgray',  # Color for Component 9
    '14': 'blue',   # Add other component numbers with specific colors
}
sc.pl.spatial(
    adata, 
    color='Component_highlight',  # Use the new column for color mapping
    palette=palette,  # Apply the custom color palette
    size=1.5,  
    spot_size=100,  
    # legend_loc='on data',  # Adjust legend location if needed
    img_key=None  # Set `None` if you don’t want the tissue image in the background
)
plt.savefig(f"/cluster/projects/schwartzgroup/vg/NEST_revision/component_{highlight_component}_plot.png", dpi=800)

