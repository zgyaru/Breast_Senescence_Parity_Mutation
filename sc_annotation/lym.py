import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import os

adata = ad.read_h5ad('../lym/adata/lym_post_scvi.h5ad')

## cell type annotation was performed in a hierarchical manner. 
## First, we annotated at a coarse level (ct_level1), code in first scvi_integration.py, 
## then we subsetted the data and performed finer annotation (ct_level2 and ct_level3)
## The code for annotating each level is shown below. 
## The final annotation is stored in ct_level3 column in the adata object.
## The annotation was based on the expression of known marker genes, the DE genes for each cluster, and the cluster dendrogram.
## dotplot of expression of known marker genes was generated using sc.pl.dotplot function, 
##and the DE genes were identified using sc.tl.rank_genes_groups function with the leiden clusters as groups.


sc.pl.umap(adata, color = ['leiden_0_3', 'leiden_0_4', 'leiden_0_5'])

sc.tl.leiden(adata, resolution=0.2, restrict_to=('leiden_0_3', ['3']), key_added='sub_0_3_3')
sc.tl.leiden(adata, resolution=0.2, restrict_to=('sub_0_3_3', ['4']), key_added='sub_0_3_3_4')
sc.tl.leiden(adata, resolution = 0.5, restrict_to=('sub_0_3_3_4', ['1']), key_added='sub_0_3_3_4_1')
sc.tl.leiden(adata, resolution = 0.3, restrict_to=('sub_0_3_3_4_1', ['6']), key_added='sub_0_3_3_4_1_6')
sc.tl.leiden(adata, resolution=0.3, restrict_to=('sub_0_3_3_4_1_6', ['2']), key_added='sub_0_3_3_4_1_6_2')
sc.tl.leiden(adata, resolution = 0.2, restrict_to=('sub_0_3_3_4_1_6_2', ['0']), key_added='sub_0_3_3_4_1_6_2_0')
sc.tl.leiden(adata, resolution = 0.2, restrict_to=('sub_0_3_3_4_1_6_2_0', ['5']), key_added='sub_0_3_3_4_1_6_2_0_5')
sc.tl.leiden(adata, resolution = 0.1, restrict_to=('sub_0_3_3_4_1_6_2_0_5', ['1,6']), key_added='sub_0_3_3_4_1_6_2_0_5_1,6')
sc.tl.dendrogram(adata, groupby='sub_0_3_3_4_1_6_2_0_5_1,6_1,3', use_rep='X_scVI', key_added='dendrogram_sub_0_3_3_4_1_6_2_0_5_1,6_1,3')
sc.tl.dendrogram(adata, groupby='sub_0_3_3_4_1_6_2_0_5_1,6_1,3_8', use_rep='X_scVI', key_added='dendrogram_sub_0_3_3_4_1_6_2_0_5_1,6_1,3_8')

cd4_dict = {
    'CD4': ['Cd4', 'Cd40lg'],
    'CD8': ['Cd8a', 'Cd8b1'],
    'CD3': ['Cd3e', 'Cd3d', 'Cd3g'],
    "CD4 Th2": ["Il3", "Il13", "Il5", "Il4", "Csf2", "Gata3"],
    "CD4 Effector": ["Ifng", "Tbx21", "Il2"],
    "CD4 Trm": ["Cd4", "Itgae"],
    "CD4 Th17": ["Il17a", "Il17f", "Il22", "Il21", "Il6"],
    "CD4 Th9": ["Irf4", "Vdr"],
    "CD4 Th22": ["Ccr10", "Il22", "Ccr4", "Ccr6"],
    "CD4 Treg": ["Foxp3"],
    "CD4 Treg activation": ["Ctla4", "Il2ra", "Gata3"],
    "CD4_Tfh": ["Bcl6", "Icos", "Pdcd1", "Cxcr5"],
}

cd8_dict = {
    'Reecptpr': ['Cd8a', 'Cd8b1', 'Cd4', 'Cd40lg'],
    'Effector': ['Gzmb', 'Klrg1'],
    'Memory': ['Il7r', 'Itgae', 'S100a4'],
    'Effector Memory': ['Il7r', 'Cx3cr1'],
    'central memory': ['Sell', 'Ccr7', 'Cd27'],
    'naive': ['Cd44'],
    'NKT': ['Il2rb','Klrb1c', 'Klra1','Klrc1','Klrc2','Klrc3','Cd3e','Cd3d','Cd3g', 'Zbtb16', 'Il10'], 
    'gammadelta': ['Rorc'],
    'Other cytotoci': ['Gzmm', 'Gzmk', 'Gzmb', 'Gzma'],
    'Exhaustion': ['Pdcd1', 'Lag3', 'Tigit', 'Havcr2'],
    'Activation': ['Cd69', 'Il2ra'],
    'memory_dis': ['Eomes', 'Tbx21'],
    'migratory': ['Cxcr3', 'Cxcr6'],
    'survive': ['Bcl2'],
    'Proliferate': ['Mki67'],
    'Effector': ['Ifng', 'Tnf']
}

level3_dict = {
    '1,3,2': 'Doublet', ##low lymphoid signature,
    '9': 'Doublet', ##high myeloid signature
    '6,3': 'Doublet', ##High Fabp4
    '1,6,1': 'Doublet', ##high fibroblast signature
    '2,5': 'Doublet', ##high stromal signature
    '7': 'Plasma', ##Doublet/plasma, high Jchain but high epi signature
    '6,5': 'Doublet', ##high myeloid signature
    '6,4': 'Doublet',#express Fabp4
    '5,2': 'Doublet', #high Fabp4
    '3,3': 'Doublet', ##low lymphoid signature, cd4 and cd8 low, brca2/palb2, exclude
    '3,4': 'Doublet', ##low lymphoid signature, cd4 and cd8 low, exclude. high mito
    '4,4': 'Doublet', ##Cd4 pos but Cd3 negative, exclude
    '6,0': 'ILC2', #cd3 cd8 cd4 negative ilc2 markers
    '6,1': 'CD4_Th2', #CD4+Il4+Il5+Th2 cells early
    '6,2': 'CD4_Th2', #CD4+Il4+Il5+Th2 cells, but lower il4 late
    '1,2': 'CD8_CTL',#CD8+Nkg7+Ccl5+ Cytotoxic T cells
    '1,9': 'Doublet', #small cluster, express both CD4 and CD8
    '1,0': 'CD8_CTL', #CD8+Ccl5+Gzmk+ Cytotoxic T cells
    '4,0': 'NK',
    '4,1': 'NK',
    '4,2': 'CD8_Tcm', #CD8+Sell+Ccr7+SellGzmm+Txk,Yes1,Ccl5Suggest early activation or functional readiness Ctla2a, Ly6c2,Mild regulatory/inflammatory profile
    '4,3': 'DN_NKT',
    '1,4': 'CD8_NKT',##express cd8a, but not cd8b1, express Nkt markers, express low levels zbtb16, but Itgae+
    '1,1': 'CD8_Trm',##cd8a+Itgae+Trm markers,
    '1,5': 'CD8_Trm/NKT',##Cd8a positive, Cd8b1 pos as well (but low), but high Gzmb, cytotoxic, only in parity likely reported by previous paper
    '1,3,3':'CD8_CTL', ##high Pdcd1, Lag3, Tigit, Havcr2, low Gzmb',
    '1,3,0':'CD8_CTL',
    '1,3,1':'CD8_CTL',
    '5,0': 'NKT17', ##high Zbtb16, Il17
    '5,1': 'NKT17', ##high Zbtb16,Il17
    '5,2': 'NKT17', ##high Zbtb16, Il17
    '10': 'CD8_Exhausted', ##high Pdcd1, Lag3, Tigit, low Gzmb',
    '2,1': 'CD4_Treg',
    '2,2': 'CD4_Effector',## CD4+Ifng+Tbx21+Th1 cells
    '2,0': 'CD4_Effector',## CD4+Ifng+Tbx21+Th1 cells
    '1,8': 'CD8_CTL',
    '2,4': 'Doublet', ##CD4 and CD8 positive,
    '1,6,0':'Proliferating_T', 
    '2,5': 'Proliferating_T',
    '0,0': 'B_cells',
    '0,1': 'B_cells',
    '0,2':'B_cells',
    '11': 'Doublet', #immature B cells remove
    '8,0': 'CD8_Tem', ##express Sell but not Ccr7 no itgae, high Cxcr3 and Cxcr6, Gzmm and Gzmk
    '8,1':'T_naive',
    '0,3': 'B_cells', 
    '1,7': 'CD8_CTL', ##high heatshock
    '2,3': 'CD4_Treg', ##express low levels of Foxp3, high heatshock 
    '3,0': 'T_naive',
    '3,1': 'T_naive',
    '3,2': 'T_naive',
    '3,5': 'Doublet', ##small cluster cd4 and cd8 negative
}


adata.obs['ct_level3'] = adata.obs['sub_0_3_3_4_1_6_2_0_5_1,6_1,3_8'].map(level3_dict)
adata.obs.loc[adata.obs['Doublet_type'] == 'Doublet', 'ct_level3'] = 'Doublet'

ct_level2 = {
    'B_cells': 'B_cells',
    'CD4_Effector': 'CD4', 
    'CD4_Th2': 'CD4',
    'CD4_Treg': 'CD4',
    'CD8_CTL': 'CD8',
    'CD8_NKT': 'CD8',
    'CD8_Trm': 'CD8',
    'CD8_Tcm': 'CD8',
    'CD8_Tem': 'CD8',
    'CD8_Exhausted': 'CD8',
    'CD8_NKT': 'NKT',
    'T_naive': 'T_naive',
    'DN_NKT': 'NKT',
    'ILC2': 'ILC',
    'NK': 'NK',
    'NKT17': 'NKT',
    'CD8_Trm/NKT': 'γδT',
    'Doublet': 'Doublet',
    'Proliferating_T': 'Proliferating_T',
    'Plasma': 'Plasma'
}

adata.obs['ct_level2'] = adata.obs['ct_level3'].map(ct_level2)

adata.write_h5ad('../lym/adata/lym_post_scvi.h5ad')



