import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import copy

# Set figure directory once
sc.settings.figdir = '../visium_HD/TACCO/wt/annotation_plot'

# Load data
wt = ad.read_h5ad('../visium_HD/Load_data/WT_adata_tacco_final.h5ad')

# Marker gene dictionary 
ct_dict = {
    'LASP_Other': ["Fos", "Atf3", "Csn2", "Csn1s1", "Csn3", "Cck", 
                   "Aldh1a3", "Ngf", "Cwh43", "Alox12e", "Krt23", "Stmn1", "Mki67"],
    'LASP7': ['Ltf', 'Ntng1'],
    'LHS': ['Pgr', 'Cited1', 'Esr1'],
    'BMYO': ['Trp63', 'Krt5', 'Krt14'],
    'T_cells': ['Cd8a', 'Cd8b1', 'Cd4', 'Cd40lg'],
    'Mac/Mono': ['Cd86', 'Adgre1'],
    'DC': ['Itgax', 'Clec9a'],
    'NK': ['Ncr1', 'Klrd1'],
    'Fibroblast': ["Col3a1", "Dcn"],
    'Vascular_cells': ["Cldn5", "Pecam1"],
    'Adipocyte': ['Fabp4', 'Adipoq'],
}

# Consolidate LASP subtypes, due to small number of spots
wt.obs['Tacco_merged'] = wt.obs['Tacco_final'].replace(
    ['LASP1', 'LASP2', 'LASP3', 'LASP4', 'LASP5', 'LASP6'], 'LASP_Other'
)

ct_order = ['LASP_Other', 'LASP7', 'LHS', 'BMYO', 'T_cells', 
            'Mac/Mono', 'DC', 'NK', 'Fibroblast', 'Vascular_cells', 'Adipocyte']

# Matrix plot
sc.pl.matrixplot(
    wt, var_names=ct_dict, groupby='Tacco_merged',
    use_raw=False, standard_scale='var',
    categories_order=ct_order,
    save='_wt_marker_genes_matrix.pdf'
)

# Palette
palette = {
    "Adipocyte": "#000000", "BMYO": "#F9F089", "DC": "#F1B883",
    "Fibroblast": "#9DCB94", "LASP7": "#CB433E", "LASP_Other": "#C27DAF",
    "LHS": "#B6A7CC", "Mac/Mono": "#7D554E", "T_cells": "#B6BC53",
    "Vascular_cells": "#9FD3DD", "NK": "#B5958F"
}

label_key = "Tacco_merged"
wt.obs[label_key] = wt.obs[label_key].astype('category')
cats = [c for c in ct_order if c in wt.obs[label_key].cat.categories]
wt.obs[label_key] = wt.obs[label_key].cat.reorder_categories(cats, ordered=True)
wt.uns[f"{label_key}_colors"] = [palette[c] for c in wt.obs[label_key].cat.categories]

# Spatial plot
fig = sc.pl.spatial(
    wt, color='Tacco_merged', alpha=1.0, alpha_img=0.0,
    palette=palette, na_color='white', size=1.5, img_key='hires',
    legend_loc='right margin', frameon=False, title='',
    show=False, return_fig=True
)
for ax in fig.axes:
    for coll in ax.collections:
        coll.set_rasterized(True)
fig.savefig(f"{sc.settings.figdir}/plot_anno.pdf", dpi=100, bbox_inches="tight")
plt.close(fig)

# Compute proportions for bar chart
counts = wt.obs[label_key].value_counts()
df = pd.DataFrame({label_key: counts.index, "proportion": counts.values / counts.sum()})
df = df.set_index(label_key).loc[cats].reset_index()  # reorder

# Stacked bar
colors = [palette[k] for k in df[label_key]]
fig, ax = plt.subplots(figsize=(3.2, 5))
bottom = 0
for ct, p, c in zip(df[label_key], df["proportion"], colors):
    ax.bar("WT", p, bottom=bottom, color=c, edgecolor="white", linewidth=0.6)
    bottom += p
ax.set_ylim(0, 1)
ax.set_ylabel("Proportion")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.legend(df[label_key], bbox_to_anchor=(1.02, 1), loc="upper left", frameon=False)
plt.tight_layout()
plt.savefig(f"{sc.settings.figdir}/wt_proportion_bar.pdf")
plt.close()

# Gene spatial plots
genes = ["Ptprc", "Krt8", "Krt14", "Cdkn2a", "Csn2", "Dcn"]
titles = ["Ptprc (CD45)", "Krt8 (K8)", "Krt14 (K14)",
          "Cdkn2a (CDKN2A)", "Csn2 (CSN2)", "Dcn (DCN)"]

white_red_cmap = copy.copy(plt.get_cmap('Reds'))
white_red_cmap.set_under('#000000')

fig = sc.pl.spatial(
    wt, color=genes, cmap=white_red_cmap, vmin=1e-2,
    show=False, img_key=None, alpha=1.0, alpha_img=0.0,
    ncols=3, size=1.5, title=titles, return_fig=True
)
for ax in fig.axes:
    for coll in ax.collections:
        coll.set_rasterized(True)
fig.savefig(f"{sc.settings.figdir}/wt_spatial_genes.pdf", dpi=100, bbox_inches="tight")
plt.close(fig)
