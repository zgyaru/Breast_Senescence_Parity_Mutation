import scanpy as sc 
import numpy as np
import os
import tacco as tc 
import pandas as pd
import anndata as ad

adata = ad.read_h5ad('../adata/adata_full_scvi.h5ad')
adata = adata[adata.obs['ct_level3'] != 'Doublet']
adata = adata[adata.obs['ct_level3'] != 'Low_Quality']
adata = adata[adata.obs['ct_level3'] != 'DDC']
adata = adata[adata.obs['ct_level3'] != 'Proliferating_T']
adata = adata[adata.obs['ct_level3'] != 'LHS1'] #both lhs and lasp marker

##remove rare cell types 


adata_anno_map={
    'LASP1': 'LASP1',
    'LASP2': 'LASP2',
    'LASP3': 'LASP3',
    'LASP4': 'LASP4',
    'LASP5': 'LASP5',
    'LASP6': 'LASP6',
    'LASP7': 'LASP7',
    'LHS2': 'LHS',
    'LHS3': 'LHS',
    'LHS4': 'LHS',
    'Fb1': 'Fibroblast',
    'Fb2': 'Fibroblast',
    'Fb3': 'Fibroblast',
    'Fb4': 'Fibroblast',
    'Fb5': 'Fibroblast',
    'Fb6': 'Fibroblast',
    'Fb7': 'Fibroblast',
    'Fb8': 'Fibroblast',
    'B_cells': 'B cell',
    'BMYO1': 'BMYO',
    'BMYO2': 'BMYO',
    'BMYO3': 'BMYO',
    'BMYO4': 'BMYO',
    'CD8_Trm': 'CD8 T cell',
    'CD8_CTL': 'CD8 T cell',
    'CD8_Tcm': 'CD8 T cell',
    'CD8_Tem': 'CD8 T cell',
    'CD8_Exhausted': 'CD8 T cell',
    'CD8_NKT': 'CD8 T cell',
    'CD4_Th1': 'CD4 T cell',
    'CD4_Th2': 'CD4 T cell',
    'CD4_Treg': 'CD4 T cell',
    'NKT17': 'NKT17',
    'NK': 'NK',
    'CD8_Trm/NKT': 'CD8 T cell',
    'Aged_Il33+': 'Aged Il33+',
    'DN_NKT': 'DN NKT',
    'Tam_1': 'Tam',
    'Tam_2': 'Tam',
    'Tam_3': 'Tam',
    'Mo2': 'Macrophage',
    'Mo1': 'Macrophage',
    'Mo3_1': 'Macrophage',
    'Mo3_2': 'Macrophage',
    'Mo3_3': 'Macrophage',
    'pDC': 'pDC',
    'mDC': 'mDC',
    'DC1': 'cDC',
    'DC2': 'cDC',
    'Non_Classical_Monocyte': 'Monocyte',
    'Classical_Monocyte': 'Monocyte',
    'Mast Cell': 'Mast Cell',
    'ILC2': 'ILC2',
    'Tumour': 'Tumour',
    'Neutrophil': 'Neutrophil',
    'T_naive': 'T naive',
    'VEA_1': 'VEA',
    'LEC_1': 'LEC',
    'VEV': 'VEV',
    'LEC_2': 'LEC',
    'PV3': 'PV',
    'PV2': 'PV',
    'VEA_2': 'VEA',
    'PV1': 'PV',
    'VEC': 'VEC',
    'VEAT': 'VEAT',
    'Plasma': 'Plasma'}

adata.obs['celltype'] = adata.obs['ct_level3'].map(adata_anno_map)

adata.obs['final_ref'] = adata.obs['celltype'].copy()
adata.obs['final_ref'] = adata.obs['final_ref'].astype(str)
##remove rare cell types and group some cell types together for better annotation of spatial data
adata = adata[~adata.obs['final_ref'].isin(['Neutrophil', 'Mast Cell', 'Aged Il33+', 'DN NKT', 'ILC2', 'T naive', 'NKT17', 'LASP6', 'LHS3', 'LHS4'])].copy()
adata.obs.loc[adata.obs['celltype'].isin(['LEC', 'VEV', 'VEA', 'VEAT', 'PV', 'VEC']), 'final_ref'] = 'Vascular_cells'
adata.obs.loc[adata.obs['celltype'].isin(['CD8 T cell', 'CD4 T cell']), 'final_ref'] = 'T_cells'
adata.obs.loc[adata.obs['celltype'].isin(['mDC', 'pDC', 'cDC']), 'final_ref'] = 'DC'
adata.obs.loc[adata.obs['celltype'].isin(['Macrophage', 'Monocyte']), 'final_ref'] = 'Mac/Mono'

adata.obs.isna().sum()

adata.write_h5ad('../visium_HD/TACCO/adata_ref.h5ad')

