import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import os

epi = ad.read_h5ad('../epi/adata/epi_scvi_post.h5ad')
epi_dict = {"Hs;Hsp": ["Pgr", "Prlr", "Esr1", "Cited1"], "LASP": ["Aldh1a3", "Fcgbp", "Lurap1l", "Csn3", "Il33", "Kit"], "LASP2": ["Elf5", "Csn2", "Wap"], "Bsl": ["Lgals7", "Trp63", "Gjb4", "Cntn2", "Krt5", "Oxtr", "Nrtn", "Arc", "Acta2", "Krt14", "Krt17"], "Tumour": ["Oas2", "Asz1"], "Stem": ['Zeb1', 'Tcf4'], 'Cycling': ['Mki67'], 'LASP5': ['Cwh43', 'Krt23', 'Cdkn2a', 'Alox12e']}
markers_dict = {"Epithelial":["Krt8","Sfn","Krt18", "Epcam"], "Fibroblast":["Col3a1","Dcn","Col4a1"],"Immune":["Ptprc","Cd52"],"Myeloid":["Tyrobp","H2-Aa"],"Lymphoid":["Hcst","Cd2","Cd3d"],"Stroma":["Pecam1", "Cdh5", "Eng", "Pdgfra", "Pdgfrb", "Fap", "Rgs5", "Des", "Notch3", "Fabp4","Apold1","Emcn","Pecam1"]}

## cell type annotation was performed in a hierarchical manner. 
## First, we annotated at a coarse level (ct_level1), code in first scvi_integration.py, 
## then we subsetted the data and performed finer annotation (ct_level2 and ct_level3)
## The code for annotating each level is shown below. 
## The final annotation is stored in ct_level3 column in the adata object.
## The annotation was based on the expression of known marker genes, the DE genes for each cluster, and the cluster dendrogram.
## dotplot of expression of known marker genes was generated using sc.pl.dotplot function, 
##and the DE genes were identified using sc.tl.rank_genes_groups function with the leiden clusters as groups.

ct_level2 = {
    '0': 'BMYO',
    '1': 'LHS',
    '2': 'LASP',
    '3': 'Doublet',
    '4': 'Tumour',
    '5': 'Aged_Il33+',
    '6': 'BMYO',
    '7': 'BMYO',
    '8': 'Doublet'
}
epi.obs['ct_level2'] = epi.obs['leiden_0_3'].map(ct_level2)
##create ct_level3 column 
epi.obs['ct_level3'] = epi.obs['ct_level2'].copy()

##BMYO

bmyo = epi[epi.obs['leiden_0_3'].isin(['6', '0', '7'])].copy()

sc.tl.leiden(bmyo, resolution = 0.5, key_added='leiden_0_5_bmyo')

level3_dict_bmyo = {
    '0': 'BMYO1',
    '1': 'BMYO1',
    '2': 'BMYO1',
    '3': 'BMYO1',
    '4': 'BMYO1',
    '5': 'BMYO2',
    '6': 'BMYO1',
    '7': 'Doublet', ##low quality
    '8': 'BMYO3',
    '9': 'BMYO4',
}

bmyo.obs['ct_level3'] = bmyo.obs['leiden_0_5_bmyo'].map(level3_dict_bmyo)

##LHS

lhs = epi[epi.obs['leiden_0_3'].isin(['1'])].copy()
sc.tl.leiden(lhs, resolution = 0.2, key_added='leiden_0_2_lhs')
sc.tl.leiden(lhs, resolution=0.2, restrict_to=('leiden_0_2_lhs', ['4']), key_added='sub_4')
sc.tl.leiden(lhs, resolution=0.3, restrict_to=('sub_4', ['0']), key_added='sub_4_0')
level3_dict_lhs = {
    '1': 'LHS1', ##lp markers
    '2':'LHS2',
    '3': 'LHS2',
    '0,1': 'LHS2',
    '0,0':'LHS2',
    '0,3': 'LHS2',
    '0,4': 'LHS2',
    '0,2': 'LHS3',
    '4,0': 'LHS4',
    '4,1': 'Doublet', ##stromal-only genes (Emb, Tcn2, Igf2/H19, Socs2).
}

lhs.obs['ct_level3'] = lhs.obs['sub_4_0'].map(level3_dict_lhs)

##LASP

lasp = epi[epi.obs['leiden_0_3'] == '2'].copy()
sc.tl.leiden(epi, resolution=0.6, key_added='leiden_0_6')
level3_dict_lasp = {
    '0': 'LASP7',
    '1': 'LASP4',
    '2': 'LASP2',
    '3': 'LASP7',
    '4': 'LASP5',
    '5': 'LASP3',
    '6': 'LASP6',
    '7': 'LASP1',
    '8': 'Low_Quality',
    '9': 'DDC', ##Sample WKAL12.3e
    '10': 'Low_Quality',
    '11': 'Doublet',
}

lasp.obs['ct_level3'] = lasp.obs['leiden_0_6'].map(level3_dict_lasp)


##Il33

l33 = epi[epi.obs['ct_level3'] == 'Aged_Il33+'].copy()

sc.tl.leiden(l33, resolution=0.3, key_added='leiden_0_3')

l33_dict = {
    '0': 'Aged_Il33+',
    '1': 'Aged_Il33+',
    '2': 'Doublet',
    '3': 'Aged_Il33+',
}
l33.obs['ct_level3'] = l33.obs['leiden_0_3'].map(l33_dict)


##extract bmyo, lhs, lasp, and l33 annotation and merge back to main adata object
bmyo_anno = bmyo.obs['ct_level3']
lhs_anno = lhs.obs['ct_level3']
lasp_anno = lasp.obs['ct_level3']
l33_anno = l33.obs['ct_level3']

epi.obs['ct_level3'] = epi.obs['ct_level3'].astype(str)
epi.obs.loc[bmyo_anno.index, 'ct_level3'] = bmyo_anno
epi.obs.loc[l33_anno.index, 'ct_level3'] = l33_anno
epi.obs.loc[lasp_anno.index, 'ct_level3'] = lasp_anno
epi.obs.loc[lhs_anno.index, 'ct_level3'] = lhs_anno

epi.write_h5ad('../epi/adata/epi_scvi_post.h5ad')