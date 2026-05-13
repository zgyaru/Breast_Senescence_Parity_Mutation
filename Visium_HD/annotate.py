import scanpy as sc
import numpy as np
import os
import tacco as tc
import pandas as pd
import anndata as ad
from pathlib import Path

# Base path for easier modification
BASE_PATH = Path('../visium_HD')

# Load data
ref_data = ad.read_h5ad(BASE_PATH / 'Load_data/BR_WT_adata_raw.h5ad')
wkbr = ad.read_h5ad(BASE_PATH / 'Load_data/WKBR_adata_raw.h5ad')
wt = ad.read_h5ad(BASE_PATH / 'Load_data/WT_adata_raw.h5ad')

# --- WKBR preprocessing ---
# Filter for barcodes > 200 counts (excludes most adipocyte regions)
wkbr = wkbr[wkbr.obs['n_counts'] > 200].copy()

# Exclude lymph node region
lymph = pd.read_csv(BASE_PATH / 'TACCO/wkbr/BR_Lymph.csv', index_col=0)
wkbr.obs = wkbr.obs.merge(lymph, left_index=True, right_index=True, how='left')
wkbr = wkbr[wkbr.obs['Lymph'] != True].copy()

# --- WT preprocessing ---
# Filter for barcodes > 200 counts
wt = wt[wt.obs['n_counts'] > 200].copy()

# Exclude region with high myocyte expression
myo = pd.read_csv(BASE_PATH / 'TACCO/wt/myo.csv', index_col=0)
wt.obs = wt.obs.merge(myo, left_index=True, right_index=True, how='left')
wt = wt[wt.obs['myo'] != 'true'].copy()

# --- Plot n_counts after filtering ---
sc.settings.figdir = str(BASE_PATH / 'TACCO/filter_barcode')
sc.pl.spatial(wkbr, color='n_counts', title='WKBR n_counts', save='_WKBR_n_counts.png')
sc.pl.spatial(wt, color='n_counts', title='WT n_counts', save='_WT_n_counts.png')

# --- Remove adipocyte clusters via k-means ---
# WKBR: remove cluster 1
k_means2_wkbr = pd.read_csv(BASE_PATH / 'TACCO/wkbr/k_means_4.csv', index_col=0)
wkbr.obs = wkbr.obs.merge(k_means2_wkbr, left_index=True, right_index=True, how='left')
wkbr = wkbr[wkbr.obs['k_means_2'] != 'Cluster 1'].copy()
sc.pl.spatial(wkbr, color='n_counts', title='WKBR n_counts', save='_WKBR_remove_kmeans1.png', color_map='coolwarm')

# WT: remove cluster 1 (FIXED: was using '4' instead of 'k_means_4')
k_means2_wt = pd.read_csv(BASE_PATH / 'TACCO/wt/k_means_2.csv', index_col=0)
wt.obs = wt.obs.merge(k_means2_wt, left_index=True, right_index=True, how='left')
wt = wt[wt.obs['k_means_2'] != 'Cluster 1'].copy()
sc.pl.spatial(wt, color='n_counts', title='WT n_counts', save='_WT_remove_kmeans1.png', color_map='coolwarm')

# --- TACCO annotation ---
# WKBR: use BRCA1 samples as reference
wkbr_ref = ref_data[ref_data.obs['Genotype'].isin(['WKBR', 'BRCA1_Tumour'])].copy()
tc.tl.annotate(wkbr, wkbr_ref, annotation_key='final_ref', result_key='tacco', counts_location='X', max_annotation=2)
wkbr.write_h5ad(BASE_PATH / 'Load_data/WKBR_adata_tacco.h5ad')

# WT: use WT nulliparous samples age < 30 as reference as WT sample is 15wks old nulliparous
wt_ref = ref_data[
    (ref_data.obs['Genotype'] == 'WT') &
    (ref_data.obs['Parity'] == 'NP') &
    (ref_data.obs['Age'] < 30)
].copy()
tc.tl.annotate(wt, wt_ref, annotation_key='final_ref', result_key='tacco', counts_location='X', max_annotation=2)
wt.write_h5ad(BASE_PATH / 'Load_data/WT_adata_tacco.h5ad')


# --- Function for extracting TACCO results ---
def get_tacco_top_annotations(df):
    """Extract top two annotations per spot from TACCO probability matrix."""
    rows = []
    for idx, row in df.iterrows():
        max_column = row.idxmax()
        max_value = row.max()
        second_max_column = row.drop(max_column).idxmax()
        second_max = row.drop(max_column).max()
        rows.append({
            'index': idx,
            'max_column': max_column,
            'max_value': max_value,
            'second_max_column': second_max_column,
            'second_max': second_max
        })
    return pd.DataFrame(rows)


# --- Extract and save TACCO results ---
# WKBR
wkbr.obsm['tacco'].to_csv(BASE_PATH / 'TACCO/wkbr/annotation_wkbr_tacco_obsm.csv')
wkbr.varm['tacco'].to_csv(BASE_PATH / 'TACCO/wkbr/markers_wkbr_tacco_varm.csv')
tacco_extracted_wkbr = get_tacco_top_annotations(wkbr.obsm['tacco'])
tacco_extracted_wkbr.to_csv(BASE_PATH / 'TACCO/wkbr/annotation_wkbr_tacco_extracted.csv')

# WT
wt.obsm['tacco'].to_csv(BASE_PATH / 'TACCO/wt/annotation_wt_tacco_obsm.csv')
wt.varm['tacco'].to_csv(BASE_PATH / 'TACCO/wt/markers_wt_tacco_varm.csv')
tacco_extracted_wt = get_tacco_top_annotations(wt.obsm['tacco'])
tacco_extracted_wt.to_csv(BASE_PATH / 'TACCO/wt/annotation_wt_tacco_extracted.csv')


# --- Finetune WKBR annotation ---
tacco_result_wkbr = pd.read_csv(BASE_PATH / 'TACCO/wkbr/annotation_wkbr_tacco_extracted.csv', index_col=0)
tacco_result_wkbr['final_anno'] = tacco_result_wkbr['max_column'].copy()

# If max is Tumour but value <= 0.8, use second best annotation
mask_low_tumour = (tacco_result_wkbr['max_value'] <= 0.8) & (tacco_result_wkbr['max_column'] == 'Tumour')
tacco_result_wkbr.loc[mask_low_tumour, 'final_anno'] = tacco_result_wkbr.loc[mask_low_tumour, 'second_max_column']

# If max is Tumour and value > 0.8, label as Tumour-like
mask_high_tumour = (tacco_result_wkbr['max_value'] > 0.8) & (tacco_result_wkbr['max_column'] == 'Tumour')
tacco_result_wkbr.loc[mask_high_tumour, 'final_anno'] = 'Tumour-like'

# Prepare final annotation dataframe
wkbr_anno_df = tacco_result_wkbr[['index', 'final_anno', 'max_value', 'second_max_column', 'second_max']].copy()
wkbr_anno_df.rename(columns={'index': 'Barcode', 'final_anno': 'Tacco'}, inplace=True)

# Mark ambiguous cells (FIXED typo: was "Umbiguous")
mask_ambiguous = (wkbr_anno_df['max_value'] - wkbr_anno_df['second_max']) <= 0.1
wkbr_anno_df['TACCO_uncertain'] = wkbr_anno_df['Tacco'].copy()
wkbr_anno_df.loc[mask_ambiguous, 'TACCO_uncertain'] = 'Ambiguous'

# Add annotation to adata
wkbr.obs = wkbr.obs.merge(
    wkbr_anno_df[['Barcode', 'Tacco', 'TACCO_uncertain']],
    left_index=True,
    right_on='Barcode',
    how='left'
)
wkbr.obs.set_index('Barcode', inplace=True)  # Restore index after merge
wkbr.write_h5ad(BASE_PATH / 'Load_data/WKBR_adata_tacco_final.h5ad')

##full path wkbr.write_h5ad('../visium_HD/Load_data/WKBR_adata_tacco_final.h5ad')


# --- Finetune WT annotation ---
tacco_result_wt = pd.read_csv(BASE_PATH / 'TACCO/wt/annotation_wt_tacco_extracted.csv', index_col=0)
tacco_result_wt['final_anno'] = tacco_result_wt['max_column'].copy()

# Mark ambiguous cells (FIXED typo: was "Umbiguous")
mask_ambiguous_wt = (tacco_result_wt['max_value'] - tacco_result_wt['second_max']) <= 0.1
tacco_result_wt['TACCO_uncertain'] = tacco_result_wt['final_anno'].copy()
tacco_result_wt.loc[mask_ambiguous_wt, 'TACCO_uncertain'] = 'Ambiguous'

wt_anno_df = tacco_result_wt[['index', 'final_anno', 'max_value', 'second_max_column', 'second_max']].copy()
wt_anno_df.rename(columns={'index': 'Barcode', 'final_anno': 'Tacco'}, inplace=True)

# Add annotation to adata
wt.obs = wt.obs.merge(
    wt_anno_df[['index', 'Tacco', 'TACCO_uncertain']],
    left_index=True,
    right_on='index',
    how='left'
)
wt.obs.drop(columns=['index'], inplace=True)
wt.write_h5ad(BASE_PATH / 'Load_data/WT_adata_tacco_final.h5ad')
