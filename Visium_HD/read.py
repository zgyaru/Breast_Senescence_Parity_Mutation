import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np
import os
import bin2cell as b2c
import anndata as ad

sc.settings.figdir = '../visium_HD/Load_data'

path = "/008um_path" ##path to the folder containing the spatial data, which should have the same structure as the output of spaceranger (with 'filtered_feature_bc_matrix' and 'spatial' subfolders)
#the image you used for --image of spaceranger, that's the one the spatial coordinates are based on
source_image_path = "../tiff_img_path"
spaceranger_image_path = "../spaceranger_img_path"
adata = b2c.read_visium(path, 
                        source_image_path = source_image_path, 
                        spaceranger_image_path = spaceranger_image_path
                       )
adata.var_names_make_unique()
adata

sc.pp.filter_genes(adata, min_cells=3)
sc.pp.filter_cells(adata, min_counts=1)

adata.layers['raw'] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.layers['log1p_norm'] = adata.X.copy()

adata.write_h5ad("../visium_HD/Load_data/hd_adata_raw.h5ad")