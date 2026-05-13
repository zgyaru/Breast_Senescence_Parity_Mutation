import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd

BASE = '../integration/final_subsets'

AGING_PATH = f'{BASE}/aging_adata.h5ad'  ##Brca1 + aging dataset
BRCA_PALB2_PATH = f'{BASE}/brca_palb2_adata.h5ad'
SENO_PATH = f'{BASE}/Seno_adata.h5ad'
WKAA_PATH = f'{BASE}/WKAA_adata.h5ad'

FULL_OUT_PATH = f'{BASE}/adata_full.h5ad'


# =========================
# Load
# =========================
aging = ad.read_h5ad(AGING_PATH)
brca_palb2 = ad.read_h5ad(BRCA_PALB2_PATH)
seno = ad.read_h5ad(SENO_PATH)
wkaa = ad.read_h5ad(WKAA_PATH)


# =========================
# Helpers
# =========================
OBS_COL_ORDER = [
    'Barcode', 'Lane_ID', 'Sample', 'Mouse_ID', 'Batch', 'Dataset',
    'Age', 'Genotype', 'Treatment', 'Doublet_score', 'Doublet_type',
    'Doublet_threshold', 'n_genes', 'n_counts', 'pct_counts_mt',
    'pct_counts_ribo', 'pct_counts_hb'
]

VAR_COL_ORDER = ['Gene_ID', 'Symbol', 'n_cells_by_counts', 'mt', 'ribo', 'hb']


def standardise_var(adata, gene_id_col, symbol_col=None):
    """
    Standardise var to:
    index = Gene_ID
    columns = Gene_ID, Symbol, n_cells_by_counts, mt, ribo, hb
    """
    adata = adata.copy()

    adata.var['Gene_ID'] = adata.var[gene_id_col]

    if symbol_col is None:
        adata.var['Symbol'] = adata.var.index.astype(str)
    else:
        adata.var['Symbol'] = adata.var[symbol_col]

    for col in VAR_COL_ORDER:
        if col not in adata.var.columns:
            adata.var[col] = np.nan

    adata.var = adata.var[VAR_COL_ORDER].copy()
    adata.var = adata.var.reset_index(drop=True)
    adata.var = adata.var.set_index('Gene_ID')
    adata.var['Gene_ID'] = adata.var.index.astype(str)
    adata.var.index.name = None

    return adata


def subset_to_gene_intersection(adatas):
    gene_sets = [set(x.var.index) for x in adatas]
    gene_intersection = set.intersection(*gene_sets)
    gene_list = sorted(gene_intersection)

    adatas_sub = [x[:, x.var.index.isin(gene_list)].copy() for x in adatas]
    return adatas_sub, gene_list


def add_raw_from_counts_if_present(adata):
    adata = adata.copy()
    if 'counts' in adata.layers:
        adata.layers['raw'] = adata.layers['counts'].copy()
    return adata


def set_obs_dtype_and_order(adata):
    adata = adata.copy()

    # safer numeric conversion
    adata.obs['Doublet_threshold'] = pd.to_numeric(adata.obs['Doublet_threshold'], errors='coerce')
    adata.obs['n_genes'] = pd.to_numeric(adata.obs['n_genes'], errors='coerce').astype('Int64')
    adata.obs['n_counts'] = pd.to_numeric(adata.obs['n_counts'], errors='coerce').astype('Int64')
    adata.obs['Age'] = pd.to_numeric(adata.obs['Age'], errors='coerce')
    adata.obs['Doublet_score'] = pd.to_numeric(adata.obs['Doublet_score'], errors='coerce')
    adata.obs['pct_counts_mt'] = pd.to_numeric(adata.obs['pct_counts_mt'], errors='coerce')
    adata.obs['pct_counts_ribo'] = pd.to_numeric(adata.obs['pct_counts_ribo'], errors='coerce')
    adata.obs['pct_counts_hb'] = pd.to_numeric(adata.obs['pct_counts_hb'], errors='coerce')

    # string columns
    for col in ['Sample', 'Mouse_ID', 'Barcode', 'Batch', 'Genotype',
                'Treatment', 'Doublet_type', 'Dataset', 'Lane_ID']:
        adata.obs[col] = adata.obs[col].astype(str)

    # ensure all expected columns exist
    for col in OBS_COL_ORDER:
        if col not in adata.obs.columns:
            adata.obs[col] = np.nan

    adata.obs = adata.obs[OBS_COL_ORDER].copy()
    return adata


def clear_unneeded_slots(adata):
    adata = adata.copy()

    # clear if present
    if hasattr(adata, "varm"):
        adata.varm.clear()
    if hasattr(adata, "obsm"):
        adata.obsm.clear()
    if hasattr(adata, "obsp"):
        adata.obsp.clear()

    return adata


# =========================
# Standardise var
# =========================
print('Standardising var')

aging = standardise_var(aging, gene_id_col='ID', symbol_col=None)
brca_palb2 = standardise_var(brca_palb2, gene_id_col='ID', symbol_col='Symbol')
seno = standardise_var(seno, gene_id_col='gene_ids', symbol_col='symbol')
wkaa = standardise_var(wkaa, gene_id_col='gene_ids', symbol_col=None)


# =========================
# Intersect genes
# =========================
print('Finding common genes')

[aging, brca_palb2, seno, wkaa], gene_list = subset_to_gene_intersection(
    [aging, brca_palb2, seno, wkaa]
)

# =========================
# Add raw layer if counts exists
# =========================
print('Adding raw layers where possible')

aging = add_raw_from_counts_if_present(aging)
brca_palb2 = add_raw_from_counts_if_present(brca_palb2)
seno = add_raw_from_counts_if_present(seno)
wkaa = add_raw_from_counts_if_present(wkaa)


# =========================
# Clean obs dtype/order
# =========================
print('Cleaning obs dtypes')

aging = set_obs_dtype_and_order(aging)
brca_palb2 = set_obs_dtype_and_order(brca_palb2)
seno = set_obs_dtype_and_order(seno)
wkaa = set_obs_dtype_and_order(wkaa)


# =========================
# Save cleaned individual objects
# =========================
aging.write_h5ad(AGING_PATH)
brca_palb2.write_h5ad(BRCA_PALB2_PATH)
seno.write_h5ad(SENO_PATH)
wkaa.write_h5ad(WKAA_PATH)


# =========================
# Concatenate
# =========================
print('Concatenating')

adatas = [aging, brca_palb2, seno, wkaa]
adata_full = ad.concat(adatas, join='inner', uns_merge = 'unique')

# remove large/unneeded slots BEFORE saving
adata_full = clear_unneeded_slots(adata_full)
adata_full.write_h5ad(FULL_OUT_PATH)

print('Done')