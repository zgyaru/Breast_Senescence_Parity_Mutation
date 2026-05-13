import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
from statannotations.Annotator import Annotator
import os


epi = ad.read_h5ad('../epi/adata/epi_scvi_post.h5ad')
stromal = ad.read_h5ad('../stromal/adata/stromal_scvi_post.h5ad')
lym = ad.read_h5ad('../lym/adata/lym_post_scvi.h5ad')
mye = ad.read_h5ad('../mye/adata/mye_scvi_post_anno.h5ad')

epi_genes = ['Cdkn2a', 'Cdkn1a', 'Cd274', 'Areg', 'Angptl4', 'Cxcl17', 'Cxcl10', 'Wnt5a'] ##LASP overall
AD_genes = ['Csn2', 'Csn1s1', 'Csn1s2a', 'Lalba', 'Wap']
# LASP5_genes = ['Cwh43', 'Alox12e', 'Krt23']
fib_genes = ['Cdkn2a', 'Cd274', 'Cxcl10'] 
lym_genes = ['Lag3', 'Tigit', 'Pdcd1'] #CD8_CTL, NK for Tigit, CTLA4 for Treg
mye_genes = ['Cd274'] ##DC1, Mo1, Mo2, Mo3_1

epi = epi[epi.obs['Group'].isin(['ABT737_control', 'Treated'])]
stromal = stromal[stromal.obs['Group'].isin(['ABT737_control', 'Treated'])]
lym = lym[lym.obs['Group'].isin(['ABT737_control', 'Treated'])]
mye = mye[mye.obs['Group'].isin(['ABT737_control', 'Treated'])]

sc.settings.figdir = '../senescence/seno'
epi = epi[~epi.obs['ct_level3'].isin(['Doublet', 'Low_Quality', 'Tumour', 'DDC'])]
LASP = epi[epi.obs['ct_level2'] == 'LASP']


epi_genes = ['Cdkn2a', 'Cdkn1a', 'Cd274', 'Areg', 'Angptl4', 'Cxcl17', 'Cxcl10', 'Wnt5a'] ##LASP overall
AD_genes = ['Csn2', 'Csn1s1', 'Csn1s2a', 'Lalba', 'Wap']
lasp_dict = {'Senescence': epi_genes, 'AD': AD_genes}

LASP_subtypes = ['LASP1', 'LASP2', 'LASP3', 'LASP4', 'LASP5', 'LASP6', 'LASP7']

for subtype in LASP_subtypes:
    LASP_sub = LASP[LASP.obs['ct_level3'] == subtype]
    sc.pl.dotplot(LASP_sub, var_names=lasp_dict, groupby='Group', use_raw=False, save=f'_ABT737_{subtype}_senad.png', show=False, color_map='Reds', standard_scale='var')
    sc.pl.dotplot(LASP_sub, var_names=lasp_dict, groupby='Group', use_raw=False, save=f'_ABT737_{subtype}_senad.svg', show=False, color_map='Reds', standard_scale='var')


dict_marker = {
    'LASP1': ["Fos", "Atf3"],
    'LASP2': ["Csn2", "Csn1s1"],
    'LASP3': ["Csn3", "Cck"],
    'LASP4': ['Aldh1a3', 'Ngf'],
    'LASP5': ['Cwh43', 'Alox12e', 'Krt23'],
    'LASP6': ['Stmn1', 'Mki67'],
    'LASP7': ['Ltf', 'Ntng1'],
}

for subtype in LASP_subtypes:
    LASP_sub = LASP[LASP.obs['ct_level3'] == subtype]
    sc.pl.dotplot(LASP_sub, var_names=dict_marker, groupby='Group', use_raw=False, save=f'_ABT737_{subtype}_marker.png', show=False, color_map='Reds', standard_scale='var')
    sc.pl.dotplot(LASP_sub, var_names=dict_marker, groupby='Group', use_raw=False, save=f'_ABT737_{subtype}_marker.svg', show=False, color_map='Reds', standard_scale='var')

sc.pl.dotplot(LASP, var_names=lasp_dict, groupby='Group', use_raw=False, save='_ABT737_LASP.png', show=False, color_map='Reds', dot_max=0.5, standard_scale='var')
sc.pl.dotplot(LASP, var_names=lasp_dict, groupby='Group', use_raw=False, save='_ABT737_epi.svg', show=False, color_map='Reds', dot_max=0.5, standard_scale='var')
sc.pl.dotplot(LASP, var_names=['Csn2', 'Csn1s1'], groupby='Group', use_raw=False, save='_ABT737_epi_csn.svg', show=False, color_map='Reds', dot_max=0.5, standard_scale='var')


fib = stromal[stromal.obs['ct_level2'] == 'Fibroblast']
fib = fib[~fib.obs['ct_level3'].isin(['Doublet'])]
fib_dict = {'Fibroblast': fib_genes}
sc.pl.dotplot(fib, var_names=fib_dict, groupby='Group', use_raw=False, save='_ABT737_fib.svg', show=False, color_map='Reds', dot_max=0.5, standard_scale='var')

lym = lym[~lym.obs['ct_level3'].isin(['Doublet'])]
CD8_CTL = lym[lym.obs['ct_level3'] == 'CD8_CTL']
NK = lym[lym.obs['ct_level3'] == 'NK']
Treg = lym[lym.obs['ct_level3'] == 'CD4_Treg']
cd8_trm = lym[lym.obs['ct_level3'] == 'CD8_Trm']

cd8_dict = {'CD8_CTL': lym_genes}
nk_dict = {'NK': lym_genes}
treg_dict = {'CD4_Treg': ['Lag3', 'Tigit', 'Pdcd1', 'Ctla4']}
cd8_trm_dict = {'CD8_Trm': lym_genes}


sc.pl.dotplot(CD8_CTL, var_names=cd8_dict, groupby='Group', use_raw=False, save='_ABT737_CD8_CTL.svg', show=False, color_map='Reds', dot_max=0.5, standard_scale='var')
sc.pl.dotplot(NK, var_names=nk_dict, groupby='Group', use_raw=False, save='_ABT737_NK.svg', show=False, color_map='Reds', dot_max=0.5, standard_scale='var')
sc.pl.dotplot(Treg, var_names=treg_dict, groupby='Group', use_raw=False, save='_ABT737_Treg.svg', show=False, color_map='Reds', dot_max=0.5, standard_scale='var')
sc.pl.dotplot(cd8_trm, var_names=cd8_trm_dict, groupby='Group', use_raw=False, save='_ABT737_CD8_Trm.svg', show=False, color_map='Reds', dot_max=0.5, standard_scale='var')


mye = mye[~mye.obs['ct_level3'].isin(['Doublet'])]
mye.obs['ct_level3'] = mye.obs['ct_level3'].astype(str)
mye.obs.loc[mye.obs['ct_level3'].isin(['Mo3_1', 'Mo3_2', 'Mo3_3']), 'ct_level3'] = 'Mo3'  # Rename Mo3 to Mo3_1 for consistency

DC1 = mye[mye.obs['ct_level3'] == 'DC1']
Mo1 = mye[mye.obs['ct_level3'] == 'Mo1']
Mo2 = mye[mye.obs['ct_level3'] == 'Mo2']
Mo3 = mye[mye.obs['ct_level3'] == 'Mo3']

dc1_dict = {'DC1': mye_genes}
mo1_dict = {'Mo1': mye_genes}
mo2_dict = {'Mo2': mye_genes}
mo3_dict = {'Mo3': mye_genes}

sc.pl.dotplot(DC1, var_names=dc1_dict, groupby='Group', use_raw=False, save='_ABT737_DC1.svg', show=False, color_map='Reds', dot_max=0.5, standard_scale='var')
sc.pl.dotplot(Mo1, var_names=mo1_dict, groupby='Group', use_raw=False, save='_ABT737_Mo1.svg', show=False, color_map='Reds', dot_max=0.5, standard_scale='var')  
sc.pl.dotplot(Mo2, var_names=mo2_dict, groupby='Group', use_raw=False, save='_ABT737_Mo2.svg', show=False, color_map='Reds', dot_max=0.5, standard_scale='var')
sc.pl.dotplot(Mo3, var_names=mo3_dict, groupby='Group', use_raw=False, save='_ABT737_Mo3_1.svg', show=False, color_map='Reds', dot_max=0.5, standard_scale='var')