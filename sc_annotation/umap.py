import scanpy as sc
import anndata as ad
import pandas as pd


epi = ad.read_h5ad('/Users/runyuxia/Desktop/HD_new_analysis/SingleCell/new_data/epi_scvi_post.h5ad')
epi = epi[epi.obs['ct_level3'] != 'Doublet']
epi = epi[epi.obs['ct_level3'] != 'Low_Quality']
epi.obs['ct_level3'] = epi.obs['ct_level3'].astype(str)


lym = ad.read_h5ad('/Users/runyuxia/Desktop/HD_new_analysis/SingleCell/new_data/lym_post_scvi.h5ad')
lym = lym[lym.obs['ct_level3'] != 'Doublet']
lym.obs['ct_level3'] = lym.obs['ct_level3'].astype(str)

mye = ad.read_h5ad('/Users/runyuxia/Desktop/HD_new_analysis/SingleCell/new_data/mye_scvi_post_anno.h5ad')
mye = mye[mye.obs['ct_level3'] != 'Doublet']
mye.obs['ct_level3'] = mye.obs['ct_level3'].astype(str)

fib = ad.read_h5ad('/Users/runyuxia/Desktop/HD_new_analysis/SingleCell/new_data/Stromal_scVI_post.h5ad')
fib = fib[fib.obs['ct_level3'] != 'Doublet']
fib.obs['ct_level3'] = fib.obs['ct_level3'].astype(str)

sc.settings.figdir = '../plot/umaps'

epi_palette ={
    'BMYO1': "#83ba26",
    "BMYO2": "#bed497",
    "BMYO3": '#b8dbca',
    "BMYO4": "#d2dc65",
    "Tumour": "#b5d7a3",
    "LASP1": '#f4abb7',
    "LASP2": "#fce1e1",
    "LASP3": '#c483b7',
    "LASP4": "#bf6688",
    "LASP5": '#9f7478',
    "LASP6": "#f18b80",
    "LASP7": "#cbbcdd",
    "LHS1": "#fcf4ad",
    'LHS2': "#ffe6a1",
    'LHS3': "#cca247",
    "LHS4": '#d3a591',
    'DDC': '#e72771',
    'Aged_Il33+': '#67b43e',
}
sc.pl.umap(epi, color='ct_level3', palette=epi_palette, frameon=False, save='_epi_ct_level3_umap.png', show=False)


lym_palette = {
    'T_naive': '#b37772',
    'CD8_Tem': '#c09462',
    'CD8_Tcm': '#bb4a4a',
    'DN_NKT': '#e2b65d',
    'NK': '#ffe5bd',
    'CD8_NKT': '#da6b4b',
    'γδT_cytotoxic': '#f6c574',
    'CD8_Trm': '#ac6e30',
    'CD8_CTL': '#e67209',
    'CD4_Effector': '#f9c561',
    'CD4_Treg': '#f4ae1f',
    'CD8_Exhausted': '#d85114',
    'CD4_Th2': '#fdeaab',
    'ILC2': '#cf9234',
    'Proliferating_T': '#8f6e67',
    'NKT17': '#b76329',
    'B_cells': '#fff4be'
}

sc.pl.umap(lym, color='ct_level3', palette=lym_palette, frameon=False, save='_lym_ct_level3_umap.png', show=False)

mye_palette = {
    'DC1': '#7b8325',
    'DC2': '#43b385',
    'Mo3_3': '#aad3a5',
    'mDC': '#9fd1ba',
    'pDC': '#a1be7e',
    'Mo1': '#a5a338',
    'Mo2': '#4db05d',
    'Classical_Monocyte': '#bed796',
    'Non_Classical_Monocyte': '#88cedd',
    'Tam_1': '#427c62',
    'Tam_2': '#6f8f61',
    'Tam_3': '#c5dfb4',
    'Mo3_1': '#5e7c31',
    'Mo3_2': '#aab685',
    'Neutrophil': '#51bdc6',
    'Mast Cell': '#dde2a1'
}
sc.pl.umap(mye, color='ct_level3', palette=mye_palette, frameon=False, save='_mye_ct_level3_umap.png', show=False)

fib_palette = {
    'Fb1': '#a24f99',
    'Fb2': '#9a7f93',
    'Fb3': '#884996',
    'Fb4': '#dec1c9',
    'Fb5': '#c9b8db',
    'Fb6': '#d293a8',
    'Fb7': '#b6a8c4',
    'Fb8': '#d9c2df',
    'LEC_1': '#bab2d9',
    'LEC_2': '#dbb5d6',
    'VEV': '#78c1cf',
    'VEAT': '#dbeaca',
    'VEC': '#6cb867',
    'VEA_1': '#7886a5',
    'VEA_2': '#b1d2f0',
    'PV1': '#cbe4cd',
    'PV2': '#828e77',
    'PV3': '#5dc5ed'
}
sc.pl.umap(fib, color='ct_level3', palette=fib_palette, frameon=False, save='_fib_ct_level3_umap.png', show=False)