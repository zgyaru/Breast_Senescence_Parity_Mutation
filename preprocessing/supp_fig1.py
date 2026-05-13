import anndata as ad
import scanpy as sc
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# ----------------------------
# Load data
# ----------------------------
adata = ad.read_h5ad(
    '../adata/adata_full_scvi.h5ad'
)

all_anno = pd.read_csv(
    '../data/all_celltype_anno.csv',
    index_col=0
)

adata.obs = adata.obs.merge(all_anno, left_index=True, right_index=True, how='left')
adata = adata[adata.obs['ct_level3'] != 'Doublet'].copy()
adata = adata[adata.obs['ct_level3'] != 'Low_Quality'].copy()
adata = adata[~adata.obs['ct_level3'].isna()].copy()

# ----------------------------
# Settings
# ----------------------------
sc.settings.figdir = '../supp_fig1/new_figs'
path = '../supp_fig1/new_figs'
output_h5ad = '../supp_fig1/adata_full_scvi_meta.h5ad'

morandi_bright_palette = [
    "#D8CAB8",
    "#BFD8B8",
    "#B8C5D8",
    "#D8B8C5",
]

print('Modifying Batch...')
batch_map = {
    '1': 'Batch 1',
    '2': 'Batch 2',
    '3': 'Batch 3',
    '4': 'Batch 4',
    '5': 'Batch 5',
    '6': 'Batch 6',
    'BRCA2_PALB2': 'Batch 7',
    'Prevention_Seno': 'Batch 8',
    'WKAA_First': 'Batch 9'
}
adata.obs['Batch_old'] = adata.obs['Batch'].copy()
adata.obs['Batch'] = adata.obs['Batch_old'].map(batch_map).fillna(adata.obs['Batch_old'])

# ----------------------------
# UMAP plots
# ----------------------------
print('Plotting UMAP before integration...')
sc.pp.neighbors(adata, n_neighbors=15, use_rep='X_pca')
sc.tl.umap(adata)
sc.pl.umap(adata, color='Batch', frameon=False, save='before_batch.png', palette='tab20')
sc.pl.umap(adata, color='Mouse_ID', frameon=False, save='before_mouse_id.png', palette='tab20')
sc.pl.umap(adata, color='Age', frameon=False, save='before_age.png', palette='tab20')
sc.pl.umap(adata, color='Parity', frameon=False, save='before_parity.png', palette='tab20')

print('Plotting UMAP after integration with scVI...')
sc.pp.neighbors(adata, n_neighbors=15, use_rep='X_scVI')
sc.tl.umap(adata)
sc.pl.umap(adata, color='Batch', palette='tab20', frameon=False, save='after_batch.png')
sc.pl.umap(adata, color='Mouse_ID', palette='tab20', frameon=False, save='after_mouse_id.png')
sc.pl.umap(adata, color='Age', palette='tab20', frameon=False, save='after_age.png')
sc.pl.umap(adata, color='Parity', palette='tab20', frameon=False, save='after_parity.png')


print('Writing adata with metadata...')
adata.write(output_h5ad)

# ----------------------------
# Reload and QC plots
# ----------------------------
adata = ad.read_h5ad(output_h5ad)

print('Plotting QC metrics...')
meta = adata.obs.copy()

meta_df = meta[['Mouse_ID', 'Batch', 'Group', 'Age', 'Parity']].drop_duplicates()
meta_df.to_csv('../supp_fig1/meta_df_ges.csv', index=False)

meta['Group'] = meta['Group'].astype('category')
meta['Group'] = meta['Group'].cat.remove_unused_categories()

meta['Mouse_ID'] = meta['Mouse_ID'].astype('category')
meta['Mouse_ID'] = meta['Mouse_ID'].cat.remove_unused_categories()

cells_per_mouse = meta.groupby('Mouse_ID', observed=True).size().reset_index(name='n_cells')
cells_per_mouse.to_csv(
    '../supp_fig1/cells_per_mouse_ges.csv',
    index=False
)

meta_batch_group = meta[['Batch', 'Mouse_ID', 'Group']].drop_duplicates()
cells_per_mouse = cells_per_mouse.merge(meta_batch_group, on='Mouse_ID', how='left')

Group_sort = [
    'Young_WT_NP', 'Old_WT_NP', 'Old_WT_Parous',
    'Young_WT_4.5dG', 'Young_WT_9.5dG', 'Young_WT_14.5dG',
    'BRCA1', 'PALB2', 'BRCA2',
    'BRCA1_Tumour', 'BRCA2_Tumour',
    'ABT737_control', 'Treated'
]

cells_per_mouse['Group'] = pd.Categorical(
    cells_per_mouse['Group'],
    categories=Group_sort,
    ordered=True
)
cells_per_mouse = cells_per_mouse.sort_values(['Group', 'Mouse_ID'])
mouse_order = cells_per_mouse['Mouse_ID'].unique().tolist()

meta['Group'] = pd.Categorical(meta['Group'], categories=Group_sort, ordered=True)
meta = meta.sort_values(['Group', 'Mouse_ID'])

palette = {
    'Young_WT_NP': '#000001',
    'Old_WT_NP': '#d28cbc',
    'Old_WT_Parous': '#6acbdd',
    'BRCA1': '#ed1f24',
    'BRCA2': '#008c45',
    'PALB2': '#f7901e',
    'BRCA1_Tumour': '#7f00ff',
    'BRCA2_Tumour': '#ff1493',
    'ABT737_control': '#1f77b4',
    'Treated': '#ffd700',
    'Young_WT_4.5dG': '#9AD9D6',
    'Young_WT_9.5dG': '#3FB7A8',
    'Young_WT_14.5dG': '#1B6F69',
}

# ----------------------------
# Plot: number of cells per mouse
# ----------------------------
print('Plotting number of cells per mouse...')
sns.set(style="whitegrid")
plt.figure(figsize=(20, 6))

ax = sns.barplot(
    data=cells_per_mouse,
    x='Mouse_ID',
    y='n_cells',
    hue='Group',
    palette=palette,
    dodge=False,
    edgecolor='black',
    width=0.8,
    order=mouse_order
)

plt.setp(ax.get_xticklabels(), rotation=90, fontsize=15)
ax.tick_params(axis='y', labelsize=15)
ax.set_xlabel('Mouse ID', fontsize=15)
ax.set_ylabel('Number of Cells', fontsize=15)
ax.set_title('Number of Cells per Mouse', fontsize=15, pad=12)
ax.legend(
    title='Group',
    bbox_to_anchor=(1.01, 1),
    loc='upper left',
    fontsize=15,
    title_fontsize=15,
    frameon=False
)
sns.despine(trim=True)
plt.tight_layout()
plt.savefig(f'{path}/cells_per_mouse_barplot.pdf', bbox_inches='tight', dpi=600)
plt.close()

# ----------------------------
# Plot: mitochondrial percentage per mouse
# ----------------------------
print('Plotting mitochondrial percentage per mouse...')
sns.set(style="whitegrid")
plt.figure(figsize=(20, 6))

meta['Mouse_ID'] = pd.Categorical(meta['Mouse_ID'], categories=mouse_order, ordered=True)

ax = sns.violinplot(
    data=meta,
    x='Mouse_ID',
    y='pct_mt_final',
    hue='Group',
    palette=palette,
    inner='box',
    scale='width',
    linewidth=1,
    cut=0,
    saturation=0.8,
    dodge=False,
    order=mouse_order
)

plt.setp(ax.get_xticklabels(), rotation=90, fontsize=15)
ax.tick_params(axis='y', labelsize=15)
ax.set_xlabel('Mouse ID', fontsize=15)
ax.set_ylabel('Mitochondrial Percentage', fontsize=15)
ax.set_title('Mitochondrial Percentage per Mouse', fontsize=15, pad=12)
ax.legend(
    title='Group',
    bbox_to_anchor=(1.01, 1),
    loc='upper left',
    fontsize=15,
    title_fontsize=15,
    frameon=False
)
sns.despine(trim=True)
plt.tight_layout()
plt.savefig(f'{path}/pct_mt_per_mouse_violinplot.pdf', bbox_inches='tight', dpi=600)
plt.close()

# ----------------------------
# Plot: number of counts per mouse
# ----------------------------
print('Plotting number of counts per mouse...')
sns.set(style="whitegrid")
plt.figure(figsize=(20, 6))

ax = sns.violinplot(
    data=meta,
    x='Mouse_ID',
    y='n_counts',
    hue='Group',
    palette=palette,
    inner='box',
    scale='width',
    linewidth=1,
    cut=0,
    saturation=0.8,
    dodge=False,
    order=mouse_order
)

ax.set_yscale('log')
plt.setp(ax.get_xticklabels(), rotation=90, fontsize=15)
ax.tick_params(axis='y', labelsize=15)
ax.set_xlabel('Mouse ID', fontsize=15)
ax.set_ylabel('Number of Counts', fontsize=15)
ax.set_title('Number of Counts per Mouse', fontsize=15, pad=12)
ax.legend(
    title='Group',
    bbox_to_anchor=(1.01, 1),
    loc='upper left',
    fontsize=15,
    title_fontsize=15,
    frameon=False
)
sns.despine(trim=True)
plt.tight_layout()
plt.savefig(f'{path}/n_counts_per_mouse_violinplot.pdf', bbox_inches='tight', dpi=600)
plt.close()

# ----------------------------
# Plot: number of genes per mouse
# ----------------------------
print('Plotting number of genes per mouse...')
sns.set(style="whitegrid")
plt.figure(figsize=(20, 6))

ax = sns.violinplot(
    data=meta,
    x='Mouse_ID',
    y='n_genes',
    hue='Group',
    palette=palette,
    inner='box',
    scale='width',
    linewidth=1,
    cut=0,
    saturation=0.8,
    dodge=False,
    order=mouse_order
)

plt.setp(ax.get_xticklabels(), rotation=90, fontsize=15)
ax.tick_params(axis='y', labelsize=15)
ax.set_ylim(bottom=0)
ax.set_xlabel('Mouse ID', fontsize=15)
ax.set_ylabel('Number of Genes', fontsize=15)
ax.set_title('Number of Genes per Mouse', fontsize=15, pad=12)
ax.legend(
    title='Group',
    bbox_to_anchor=(1.01, 1),
    loc='upper left',
    fontsize=15,
    title_fontsize=15,
    frameon=False
)
sns.despine(trim=True)
plt.tight_layout()
plt.savefig(f'{path}/n_genes_per_mouse_violinplot.pdf', bbox_inches='tight', dpi=600)
plt.close()
