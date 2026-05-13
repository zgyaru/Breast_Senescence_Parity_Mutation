import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import os

stro = ad.read_h5ad('../stromal/adata/stromal_scvi_post.h5ad')

## cell type annotation was performed in a hierarchical manner. 
## First, we annotated at a coarse level (ct_level1), code in first scvi_integration.py, 
## then we subsetted the data and performed finer annotation (ct_level2 and ct_level3)
## The code for annotating each level is shown below. 
## The final annotation is stored in ct_level3 column in the adata object.
## The annotation was based on the expression of known marker genes, the DE genes for each cluster, and the cluster dendrogram.
## dotplot of expression of known marker genes was generated using sc.pl.dotplot function, 
## and the DE genes were identified using sc.tl.rank_genes_groups function with the leiden clusters as groups.

sc.tl.leiden(stro, resolution = 0.2, restrict_to=('leiden_0_4', ['4']), key_added='sub_4')
sc.tl.leiden(stro, resolution = 0.2, restrict_to=('sub_4', ['5']), key_added='sub_4_5')

stro_dict = {
    "Vascular mural" : ["Notch3", "Myh11", "Pdgfrb"],
    "VM1": ["Adra2c", "Ctsc"],
    "VM2": ["Gpat2", "Pgap1"],
    "VM3": ["Rergl", "Acta2"],
    "VM4": ["Adgrl3", "Ggt5"], 
    "VM5": ["Rgs5", "Thbs4"],
    "EC": ["Plvap", "Cldn5", "Pecam1"],
    "EC venous": ["Ackr1", "Selp"],
    "EC capillary": ["Cd36", 'Plvap', 'Emcn'],
    "EC arterial": ["Sox17", "Efnb2", "Dll4", "Efnb2"],
    "EC angiogenic tip": ["Pxdn", "Angpt2", "Esm1"],
    "Lymphatic EC": ["Mmrn1", "Pdpn"],
    "Lymphatic EC1": ["Sema3d", "B3gnt7"],
    "Lymphatic EC2": ["Lyve1", "Tbx1"],
    "Epithelial":["Krt8","Sfn","Krt18", "Epcam"], "Fibroblast":["Col3a1","Dcn","Col4a1"],"Immune":["Ptprc","Cd52"],"stroloid":["Tyrobp","H2-Aa"],"Lymphoid":["Hcst","Cd2","Cd3d"],"Stroma":["Pecam1", "Cdh5", "Eng", "Pdgfra", "Pdgfrb", "Fap", "Rgs5", "Des", "Notch3", "Fabp4","Apold1","Emcn","Pecam1"],
    "Swc": ['Plp1','Mbp', 'Mpz'],
    "strosenchymal": ['Cd44', 'Nt5e', 'Pdgfra', 'Itgb1', 'Eng', 'Vcam1'],
    'Pan': ['Pecam1', 'Kdr', 'Cldn5', 'Emcn', 'Cdh5', 'Tie1', 'Egfl7'],
    'Vein': ['Nr2f2', 'Vwf'],
    'Aetery': ['Gja4', 'Mecom'],
    "Epithelial":["Krt8","Sfn","Krt18"], 
    "Fibroblast":["Col3a1","Dcn","Col4a1"],
    "Immune":["Ptprc","Cd52"],
    "stroloid":["Tyrobp","H2-Aa"],
    "Lymphoid":["Hcst","Cd2","Cd3d"],
    "Stroma":["Fabp4","Apold1","Emcn","Pecam1"]
}

##FIB
fib = stro[stro.obs['leiden_0_3'].isin(['1', '7', '6', '10', '0', '2', '4', '9', '12', '11'])]
sc.tl.leiden(fib, resolution = 0.2, key_added='fib_leiden_0_2')
sc.tl.leiden(fib, resolution = 0.2, restrict_to=('fib_leiden_0_2', ['4']), key_added='sub_0_2_4')

fib_dict = {
    '5': 'Doublet',
    '8': 'Fb7', ##Hmcn2	ECM glycoprotein; associated with basement membrane zones
    '6': 'Fb4',
    '12': 'Doublet',
    '1': 'Fb1', ##Pi16 pos
    '2': 'Fb2', ##Pi16 neg, Lpl, Fabp4 pos
    '3': 'Fb3', ##Gdf10 pos
    '4,2': 'Fb8',##proliferating
    '4,0': 'Fb5',
    '4,1': 'Fb5',
    '4,3': 'Fb5',
    '4,4': 'Fb5', 
    '4,5': 'Doublet', ##low n_counts, low quality
    '0': 'Fb5', ##cluster 4 and 0 pretty similar Enpp2, Slit2 pos,
    '7': 'Fb6', ##Ccl11 and Crispld2 suggest these fibroblasts may also participate in immune signaling, pointing toward inflammatory or activated fibroblasts (e.g., iCAFs).Apod and Cp might suggest stress response or senescence-associated secretory phenotype (SASP) components.
    '9': 'Doublet', ##low n_counts, low quality
    '10': 'Doublet',##Pi16 pos, Cdh5 pos, doublet
    '11': 'Doublet' ##low n_counts, low quality
}

fib.obs['ct_level2'] = 'Fibroblast'
fib.obs['ct_level3'] = fib.obs['sub_0_2_4'].map(fib_dict)
fib.obs.loc[fib.obs['Doublet_type'] == 'Doublet', 'ct_level3'] = 'Doublet'

##VASCULAR

vas = stro[~stro.obs['leiden_0_3'].isin(['1', '7', '6', '10', '0', '2', '4', '9', '12', '11'])]
sc.tl.leiden(vas, resolution = 0.4, key_added='vas_leiden_0_4')
sc.tl.leiden(vas, resolution = 0.2, restrict_to=('vas_leiden_0_4', ['6']), key_added='vas_sub_6')
sc.tl.leiden(vas, resolution = 0.2, restrict_to=('vas_sub_6', ['2']), key_added='vas_sub_6_2')
sc.tl.leiden(vas, resolution = 0.2, restrict_to=('vas_sub_6_2', ['10']), key_added='vas_sub_6_2_10')
sc.tl.leiden(vas, resolution = 0.2, restrict_to=('vas_sub_6_2_10', ['4']), key_added='vas_sub_6_2_10_4')

vas_dict = {
    '11': 'Doublet', ##Dcn
    '12': 'Doublet', 
    '8': 'Doublet',
    '4': 'Doublet',
    '10,0': 'VEAT',
    '10,1': 'Doublet',
    '10,2': 'VEV',
    '10,3':'VEAT',
    '6,2': 'LEC_2', ##has the shape of a doublet lyve1 pos
    '0': 'VEA_1',
    '1': 'VEV',
    '2,0': 'PV1',
    '2,1': 'PV2',
    '5': 'PV3',
    '3': 'VEA_2',
    '6,0': 'LEC_1',
    '6,1': 'LEC_2',
    '7': 'VEC',
    '9': 'Doublet'   
}

vas.obs['ct_level2'] = 'Vascular'
vas.obs['ct_level3'] = vas.obs['vas_sub_6_2_10'].map(vas_dict)
vas.obs.loc[vas.obs['Doublet_type'] == 'Doublet', 'ct_level3'] = 'Doublet'


##extract stromal cells and combine annotation
vas_anno = vas.obs['ct_level3']
fib_anno = fib.obs['ct_level3']

stro.obs['ct_level2'] = 'Stroma'
stro.obs['ct_level3'] = 'Stroma'

stro.obs['ct_level3'] = stro.obs['ct_level3'].astype(str)
stro.obs.loc[vas_anno.index, 'ct_level3'] = vas_anno
stro.obs.loc[fib_anno.index, 'ct_level3'] = fib_anno

stro.obs.loc[vas_anno.index, 'ct_level2'] = 'Vascular'
stro.obs.loc[fib_anno.index, 'ct_level3'] = 'Fibroblast'

stro.write_h5ad('../stro/adata/stro_scvi_post.h5ad')


