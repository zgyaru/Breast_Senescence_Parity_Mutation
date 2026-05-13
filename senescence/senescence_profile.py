import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

epi = ad.read_h5ad('../epi/adata/epi_scvi_post.h5ad')
epi.obs['Group'] = epi.obs['Group'].astype(str)
epi = epi[epi.obs['ct_level2'].isin(['LASP'])]
epi = epi[epi.obs['ct_level3'] != 'Doublet']
epi = epi[epi.obs['ct_level3'] != 'DDC']
epi = epi[epi.obs['ct_level3'] != 'Low_Quality']
epi = epi[epi.obs['Parity'].isin(['NP', 'Parous'])]
epi = epi[~epi.obs['Group'].isin(['ABT737_control', 'Treated'])]

marker_dict = {
    "LASP2": [
        "Isg15", "Ifit1", "Ifit2", "Ifit3", "Oas1g",
        "Stat1", "Stat2", "Irf9", "Il6", "Il1a", "Il1b"
    ],
    "LASP5": [
        "Tnf", "Ccl2", "Cxcl1", "Cxcl2", "Tnfaip3", "Nfkbia", "Nfkbiz"
    ]
}

sc.pl.dotplot(epi, var_names=marker_dict, groupby='ct_level3', standard_scale='var', save='_lasp2_lasp5_sen_dotplot_level3.png', use_raw=False, dot_max = 0.5)