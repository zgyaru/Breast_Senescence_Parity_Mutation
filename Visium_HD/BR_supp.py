import scanpy as sc
import pandas as pd
import anndata as ad
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

adata = ad.read_h5ad('../visium_HD/Load_data/WKBR_adata_tacco_final.h5ad')

sc.settings.figdir = '../visium_HD/TACCO/wkbr/annotation_plot'

ct_dict = {
    'Epithelial': ["Krt8", "Krt18", "Epcam"],
    'LASP': ['Kit', 'Fcgbp'],
    'LASP1': ["Fos", "Atf3"],
    'LASP2': ["Csn2", "Csn1s1"],
    'LASP3': ["Csn3", "Cck"],
    'LASP4': ['Aldh1a3', 'Ngf'],
    'LASP5': ['Cwh43', 'Alox12e', 'Krt23'],
    'LASP6': ['Stmn1', 'Mki67'],
    'LASP7': ['Ltf', 'Ntng1'],
    'LHS': ['Pgr', 'Cited1', 'Esr1'],
    'BMYO': ['Trp63', 'Krt5', 'Krt14'],
    'Tumour': ['Tk1', 'Stmn1'],
    'Immune': ['Ptprc', 'Cd52'],
    'CD8': ['Cd8a', 'Cd8b1'],
    'CD4': ['Cd4', 'Cd40lg'],
    'B_cells': ['Cd19', 'Cd79a'],
    'Plasma': ['Jchain'],
    'Macrophage': ['Cd86', 'Adgre1'],
    'Tam': ['Spp1'],
    'Dendritic': ['Itgax', 'Clec9a'],
    'Fibroblast': ["Col3a1", "Dcn"],
    "Vascular": ["Cldn5", "Pecam1"],
    'Adipocyte': ['Fabp4', 'Adipoq']
}

# Check actual column name - adjust if needed
# print(adata.obs.columns.tolist())  # Uncomment to verify

# Merge rare LASP subtypes (use correct source column name)
adata.obs['Tacco_merged'] = adata.obs['Tacco'].replace(
    ['LASP1', 'LASP3', 'LASP4', 'LASP6', 'LASP7'], 'LASP_Other'
)

# Also merge for uncertain annotations if needed
adata.obs['TACCO_uncertain_merged'] = adata.obs['TACCO_uncertain'].replace(
    ['LASP1', 'LASP3', 'LASP4', 'LASP6', 'LASP7'], 'LASP_Other'
)

# Category order - must match actual values in the column
# Verify with: print(adata.obs['Tacco_merged'].unique())
ct_order = [
    'LASP2', 'LASP5', 'LASP_Other', 'LHS', 'BMYO', 'Tumour-like',
    'CD8', 'CD4', 'B_cells', 'Plasma', 'Macrophage', 'Tam', 'Dendritic',
    'Fibroblast', 'Vascular', 'Adipocyte'
]

# Filter ct_order to only include categories that exist in the data
existing_categories = adata.obs['Tacco_merged'].unique()
ct_order_filtered = [c for c in ct_order if c in existing_categories]

# Plots
sc.pl.matrixplot(
    adata, var_names=ct_dict, groupby='Tacco_merged', use_raw=False,
    standard_scale='var', categories_order=ct_order_filtered,
    save='_wkbr_marker_genes_matrix.png'
)

sc.pl.matrixplot(
    adata, var_names=ct_dict, groupby='Tacco_merged', use_raw=False,
    standard_scale='var', categories_order=ct_order_filtered,
    save='_wkbr_marker_genes_matrix.pdf'
)

sc.pl.matrixplot(
    adata, var_names=ct_dict, groupby='Tacco_merged', use_raw=False,
    standard_scale='var', categories_order=ct_order_filtered,
    save='_wkbr_marker_genes_matrix.svg'
)
