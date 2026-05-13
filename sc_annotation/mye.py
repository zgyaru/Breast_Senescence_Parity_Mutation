import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import os

mye = ad.read_h5ad('../mye/adata/mye_scvi_post.h5ad')
mye.X = mye.layers['log1p_norm']

## cell type annotation was performed in a hierarchical manner. 
## First, we annotated at a coarse level (ct_level1), code in first scvi_integration.py, 
## then we subsetted the data and performed finer annotation (ct_level2 and ct_level3)
## The code for annotating each level is shown below. 
## The final annotation is stored in ct_level3 column in the adata object.
## The annotation was based on the expression of known marker genes, the DE genes for each cluster, and the cluster dendrogram.
## dotplot of expression of known marker genes was generated using sc.pl.dotplot function, 
##and the DE genes were identified using sc.tl.rank_genes_groups function with the leiden clusters as groups.

mye_dict = mye_dict = {'Myeloid':["Tyrobp","H2-Aa"],'DC': ['Fcgr1', 'Adgre1', 'Mertk', 'H2-Ab1', 'Itgax', 'Cd80', 'Cd86', 'Cd40', 'H2-Eb1', 'H2-Eb2', 'Flt3', 'Zbtb46'],'DC1': ['Irf8', 'Itgae'], 'DC2': ['Itgam', 'Esam', 'Tbx21', 'Rorc'], 'Macrophage': ['Fcgr1', 'Adgre1', 'Mertk', 'Csf1'], 'Mo1/Mo2': ['Mrc1', 'Itgam'], 'Mo1':['Lyve1', 'Cd86'], 'Mo2': ['Lyve1', 'Il10', 'Pdcd1lg2', 'Il4'], 'Mo3': ['Itgax', 'Cx3cr1', 'Tmem119', 'Lyve1', 'Hexb', 'Lgmn', 'Mmp12', 'Mmp13', 'Mmp14', 'Ctsf'], 'Tam': ['Itgax', 'Vcam1', 'Vegfa'], 'Mast cell': ['Cpa3', 'Tpsb2'],
 'Neutrophil': ['S100a9', 'S100a8'], 'mDC': ['Fscn1', 'Ccl22', 'Il4i1'],'pDC ': ['Ptgds', 'Siglech', 'Irf8', 'Tcf4'], 'Eosinphil': ['Siglecf', 'Il5ra'], 'Basophil': ['Il3ra', 'H2-Ab1'], 'B_cell': ['Cd19'], 'Monocyte': ['Ly6c1', 'Ly6c2', 'Itgam', 'Csf1r'], 'Plasma' : ['Il6ra', 'Sdc1'], 'Adaptive': ['Ptprc'], 'mregDC': ['Ccr7', 'Fscn1', 'Cd274'], 'M1': ['Cd86', 'Cd38', 'Il1b', 'Il6', 'Tnf'], 'M2':['Cd163', 'Mrc1'], 'cDC1': ['Clec9a', 'Xcr1', 'Cd8a'], 'cDC2': ['Clec10a', 'Fcer1a', 'Irf4', 'Cd4'], 'myephoid': ["Hcst","Cd2","Cd3d"], "Epithelial":["Krt8","Sfn","Krt18"], "Fibroblast":["Col3a1","Dcn","Col4a1"]}
markers_dict = {"Epithelial":["Krt8","Sfn","Krt18", "Epcam"], "Fibroblast":["Col3a1","Dcn","Col4a1"],"Immune":["Ptprc","Cd52"],"Myeloid":["Tyrobp","H2-Aa"],"Lymphoid":["Hcst","Cd2","Cd3d"],"Stroma":["Pecam1", "Cdh5", "Eng", "Pdgfra", "Pdgfrb", "Fap", "Rgs5", "Des", "Notch3", "Fabp4","Apold1","Emcn","Pecam1"], "B_cell": ['Blnk', 'Cd79a', 'Cd79b'], "Plasma": ['Jchain']}

sc.tl.leiden(mye, resolution=0.2, restrict_to=('leiden_0_6', ['6']), key_added='sub_6')
sc.tl.leiden(mye, resolution=0.5, restrict_to=('sub_6', ['0']), key_added='sub_6_0')
sc.tl.leiden(mye, resolution=0.2, restrict_to=('sub_6_0', ['2']), key_added='sub_6_0_2')
sc.tl.leiden(mye, resolution=0.2, restrict_to=('sub_6_0_2', ['6,2']), key_added='sub_6_0_2_6,2')

mye_dict_level2 = {
    '9': 'Doublet',
    '6,3': 'Doublet',
    '13': 'Doublet',
    '15': 'Doublet',
    '6,2,1': 'Doublet',
    '0,0': 'Macrophage',
    '0,1': 'Macrophage',
    '0,2': 'Macrophage',
    '0,3': 'Macrophage',
    '0,4': 'Macrophage',
    '0,5': 'Macrophage',
    '1': 'Macrophage',
    '2,0': 'Monocyte',
    '2,1': 'Monocyte',
    '3': 'DC',
    '4': 'Neutrophil',
    '5': 'Tam',
    '6,0': 'Macrophage',
    '6,1': 'Macrophage',
    '6,2,0': 'Tam',
    '7': 'Macrophage',
    '8': 'Tam',
    '10': 'DC',
    '11': 'DC',
    '12': 'Mast Cell',
    '14': 'DC'
}

mye.obs['ct_level2'] = mye.obs['sub_6_0_2_6,2'].map(mye_dict_level2)
mye.obs.loc[mye.obs['Doublet_type'] == 'Doublet', 'ct_level2'] = 'Doublet'

mye_dict_level3 = {
    '9': 'Doublet',
    '6,3': 'Doublet',
    '13': 'Doublet',
    '15': 'Doublet',
    '6,2,1': 'Doublet',
    '0,0': 'Mo2',
    '0,1': 'Mo1',
    '0,2': 'Mo2',
    '0,3': 'Mo2',
    '0,4': 'Mo2',
    '0,5': 'Mo2',
    '1': 'Mo3_1',
    '2,0': 'Classical_Monocyte',
    '2,1': 'Non_Classical_Monocyte',
    '3': 'DC2',
    '4': 'Neutrophil',
    '5': 'Tam_1',
    '6,0': 'Mo3_2',
    '6,1': 'Mo3_2',
    '6,2,0': 'Tam_2',
    '7': 'Mo3_3',
    '8': 'Tam_3',
    '10': 'mDC',
    '11': 'DC1',
    '12': 'Mast Cell',
    '14': 'pDC'
}

mye.obs['ct_level3'] = mye.obs['sub_6_0_2_6,2'].map(mye_dict_level3)
mye.obs.loc[mye.obs['Doublet_type'] == 'Doublet', 'ct_level3'] = 'Doublet'

mye.write_h5ad('../mye/adata/mye_scvi_post.h5ad')