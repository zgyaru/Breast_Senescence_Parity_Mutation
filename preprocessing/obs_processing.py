import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd


# =========================
# Paths
# =========================
PREVENTION_PATH = '../adata_with_doublet_scores.h5ad'
AGING_PATH = '../adata_with_doublet_scores_aging.h5ad'
WKAA_PATH = '../adata_with_doublet_scores_wkaa.h5ad'
BRCA_PALB2_PATH = '../adata_with_doublet_scores_palb2.h5ad'
AGING_META_PATH = '../data_subset_meta/aging_meta.csv'

OUTDIR = '../final_subsets'
FIGDIR = '../QC_results'

sc.settings.figdir = FIGDIR


# =========================
# Helpers
# =========================
FINAL_OBS_COLS = [
    'Parity', 'Doublet_threshold', 'Genotype_full', 'Sample', 'Mouse_ID',
    'Barcode', 'Batch', 'Age', 'Genotype', 'Treatment', 'Doublet_score',
    'Doublet_type', 'n_counts', 'n_genes', 'pct_counts_mt', 'pct_counts_ribo',
    'Dataset', 'Lane_ID', 'pct_counts_hb'
]


def ensure_obs_columns(adata, columns, fill_value=np.nan):
    """Make sure all requested obs columns exist before subsetting."""
    for col in columns:
        if col not in adata.obs.columns:
            adata.obs[col] = fill_value
    return adata


def standardise_seno_age(age_series):
    """
    Convert strings like 46W3D / 46W into week numbers.
    If days > 4, round up by 1 week.
    """
    tmp = age_series.astype(str).str.extract(r'(?P<Weeks>\d+)W(?:(?P<Days>\d+)D)?')
    tmp['Weeks'] = pd.to_numeric(tmp['Weeks'], errors='coerce')
    tmp['Days'] = pd.to_numeric(tmp['Days'], errors='coerce').fillna(0)
    return tmp['Weeks'] + (tmp['Days'] > 4).astype(int)


def calculate_qc_if_needed(adata):
    """Calculate mt/ribo/hb QC metrics if pct columns are missing."""
    needed = ['pct_counts_mt', 'pct_counts_ribo', 'pct_counts_hb']
    if all(col in adata.obs.columns for col in needed):
        return adata

    adata = adata.copy()
    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl"))
    adata.var["hb"] = adata.var_names.str.contains("^Hb[^(p)]", regex=True)

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo", "hb"],
        inplace=True,
        log1p=True
    )
    return adata


def safe_order_group_categories(adata, order):
    existing = adata.obs['Group'].dropna().astype(str).unique().tolist()
    categories = [x for x in order if x in existing] + [x for x in existing if x not in order]
    adata.obs['Group'] = pd.Categorical(adata.obs['Group'].astype(str), categories=categories, ordered=True)
    return adata


# =========================
# Load data
# =========================
prevention_adata = ad.read_h5ad(PREVENTION_PATH)
aging_adata = ad.read_h5ad(AGING_PATH)
WKAA_adata = ad.read_h5ad(WKAA_PATH)
brca_palb2 = ad.read_h5ad(BRCA_PALB2_PATH)


# =========================
# 1. Clean Senolytics
# =========================
Seno = prevention_adata[
    prevention_adata.obs['test_group'].isin(['Seno_ABT737_control', 'Seno_ABT737'])
].copy()

Seno.obs['Dataset'] = 'Senolytics'
Seno.obs['Genotype_full'] = Seno.obs['genotype']
Seno.obs['Treatment'] = Seno.obs['test_group']
Seno.obs['Age'] = Seno.obs['age']
Seno.obs['Batch'] = 'Prevention_Seno'
Seno.obs['Lane_ID'] = Seno.obs['lane_id']
Seno.obs['Barcode'] = Seno.obs['barcode']
Seno.obs['Doublet_score'] = Seno.obs['doublet_score']
Seno.obs['Doublet_type'] = Seno.obs['doublet_status']
Seno.obs['Doublet_threshold'] = Seno.obs['doublet_threshold']
Seno.obs['Sample'] = Seno.obs['sample_name_cellranger']
Seno.obs['Mouse_ID'] = Seno.obs['mouse_id']
Seno.obs['Genotype'] = 'WKAA'
Seno.obs['Parity'] = 'NP'
Seno.obs['Condition'] = 'Senolytics'

Seno.obs['Age'] = standardise_seno_age(Seno.obs['Age'])

Seno = ensure_obs_columns(Seno, FINAL_OBS_COLS)
Seno.obs = Seno.obs[FINAL_OBS_COLS]

sc.pl.umap(Seno, color='Treatment', save='Seno_umap_test_group.png')
Seno.write_h5ad(f'{OUTDIR}/Seno_adata.h5ad')


# =========================
# 2. Clean WKAA
# =========================
WKAA_adata = calculate_qc_if_needed(WKAA_adata)

WKAA_adata.obs['Dataset'] = 'WKAA_First'
WKAA_adata.obs['Genotype_full'] = 'BLG-Cre:Hemi; Brca1:Hom; p53flox:Het'
WKAA_adata.obs['n_counts'] = WKAA_adata.obs['total_counts']
WKAA_adata.obs['n_genes'] = WKAA_adata.obs['n_genes_by_counts']
WKAA_adata.obs['Treatment'] = 'None'
WKAA_adata.obs['Batch'] = 'WKAA_First'
WKAA_adata.obs['Lane_ID'] = 'None'
WKAA_adata.obs['Barcode'] = WKAA_adata.obs.index.astype(str)
WKAA_adata.obs['Doublet_score'] = WKAA_adata.obs['doublet_score']
WKAA_adata.obs['Doublet_type'] = WKAA_adata.obs['doublet_status']
WKAA_adata.obs['Doublet_threshold'] = WKAA_adata.obs['doublet_threshold']
WKAA_adata.obs['Mouse_ID'] = WKAA_adata.obs['mouse_ID']
WKAA_adata.obs['Genotype'] = 'WKAA'
WKAA_adata.obs['Parity'] = 'NP'

# Fill if missing
if 'Sample' not in WKAA_adata.obs.columns:
    WKAA_adata.obs['Sample'] = WKAA_adata.obs['Mouse_ID']
if 'Age' not in WKAA_adata.obs.columns:
    WKAA_adata.obs['Age'] = np.nan

WKAA_adata = ensure_obs_columns(WKAA_adata, FINAL_OBS_COLS)
WKAA_adata.obs = WKAA_adata.obs[FINAL_OBS_COLS]

sc.pl.umap(WKAA_adata, color='Sample', save='WKAA_umap_sample.png')
WKAA_adata.write_h5ad(f'{OUTDIR}/WKAA_adata.h5ad')


# =========================
# 3. Clean Aging atlas
# =========================
aging_adata = aging_adata[aging_adata.obs['Batch'].isin(['1', '2', '3', '4', '5', '6'])].copy()
aging_adata = calculate_qc_if_needed(aging_adata)

aging_adata.obs['Parity'] = 'NP'
aging_adata.obs.loc[aging_adata.obs['Condition'] == 'ParOld', 'Parity'] = 'Parous'
aging_adata.obs.loc[aging_adata.obs['Condition'].isin(['4.5dG', '9.5dG', '14.5dG']), 'Parity'] = 'G'

aging_adata.obs['Genotype'] = 'WT'
aging_adata.obs.loc[
    aging_adata.obs['Group'].isin(['WKBR_Tumour', 'WKBR_MG']),
    'Genotype'
] = 'WKBR'

aging_adata.obs['Dataset'] = 'Aging_atlas'
aging_adata.obs['n_counts'] = aging_adata.obs['total_counts']
aging_adata.obs['n_genes'] = aging_adata.obs['n_genes_by_counts']
aging_adata.obs['Mouse_ID'] = aging_adata.obs['sample.IDs']
aging_adata.obs['Treatment'] = 'None'
aging_adata.obs['Lane_ID'] = 'None'
aging_adata.obs['Barcode'] = aging_adata.obs.index.astype(str)
aging_adata.obs['Doublet_score'] = aging_adata.obs['doublet_score']
aging_adata.obs['Doublet_type'] = aging_adata.obs['doublet_status']
aging_adata.obs['Doublet_threshold'] = aging_adata.obs['doublet_threshold']
aging_adata.obs['Genotype_full'] = aging_adata.obs['Genotype']

# Fill if missing
if 'Sample' not in aging_adata.obs.columns:
    aging_adata.obs['Sample'] = aging_adata.obs['Mouse_ID']
if 'Age' in aging_adata.obs.columns:
    aging_adata.obs['Age'] = pd.to_numeric(aging_adata.obs['Age'], errors='coerce')
else:
    aging_adata.obs['Age'] = np.nan

aging_adata = ensure_obs_columns(aging_adata, FINAL_OBS_COLS)
aging_adata.obs = aging_adata.obs[FINAL_OBS_COLS]

sc.pl.umap(aging_adata, color='Batch', save='aging_umap_batch.png')
aging_adata.write_h5ad(f'{OUTDIR}/aging_adata.h5ad')


# =========================
# 4. Clean BRCA2 / PALB2
# =========================
brca_palb2 = calculate_qc_if_needed(brca_palb2)

brca_palb2.obs['Dataset'] = 'BRCA2_PALB2'
brca_palb2.obs['Genotype_full'] = brca_palb2.obs['KO']
brca_palb2.obs['Mouse_ID'] = brca_palb2.obs['Sample']
brca_palb2.obs['Treatment'] = 'None'
brca_palb2.obs['Batch'] = 'BRCA2_PALB2'
brca_palb2.obs['Lane_ID'] = 'None'
brca_palb2.obs['Barcode'] = brca_palb2.obs['Barcode']
brca_palb2.obs['Doublet_score'] = brca_palb2.obs['doublet_score']
brca_palb2.obs['Doublet_type'] = brca_palb2.obs['doublet_status']
brca_palb2.obs['Doublet_threshold'] = brca_palb2.obs['doublet_threshold']
brca_palb2.obs['Genotype'] = brca_palb2.obs['KO']
brca_palb2.obs['Parity'] = 'NP'

# Add BRCA/PALB2 ages
age_df = pd.DataFrame({
    'Sample': [
        'WKAD11.2e', 'WKAD11.2f', 'WKAD11.3e', 'WKAD17.1g', 'WKAD18.2f',
        'WKAL12.1e', 'WKAL12.3e', 'WKAL12.4j', 'WKAL14.1c', 'WKAL8.2j',
        'WKAL8.3l', 'WKWT41.2f'
    ],
    'Age': [46.1, 46.1, 43, 32.1, 26.4, 44.3, 37.1, 34, 24.2, 47.3, 41.2, 22.1]
})

brca_palb2.obs['Age'] = brca_palb2.obs['Sample'].map(age_df.set_index('Sample')['Age'])
brca_palb2.obs['Age'] = pd.to_numeric(brca_palb2.obs['Age'], errors='coerce')

brca_palb2 = ensure_obs_columns(brca_palb2, FINAL_OBS_COLS)
brca_palb2.obs = brca_palb2.obs[FINAL_OBS_COLS]

sc.pl.umap(brca_palb2, color='Genotype', save='brca_palb2_umap_KO.png')
brca_palb2.write_h5ad(f'{OUTDIR}/brca_palb2_adata.h5ad')


# =========================
# 5. Build parity / tumour reference from aging metadata
# =========================
aging_meta = pd.read_csv(AGING_META_PATH)

aging_meta['Parity'] = 'NP'
aging_meta.loc[aging_meta['Condition'] == 'ParOld', 'Parity'] = 'Parous'
aging_meta.loc[aging_meta['Condition'].isin(['4.5dG', '9.5dG', '14.5dG']), 'Parity'] = 'G'

parity_info = aging_meta[['Sample', 'Parity', 'Group']].drop_duplicates()

tumour_samples = parity_info.loc[parity_info['Group'] == 'WKBR_Tumour', 'Sample'].unique().tolist()
parous_samples = parity_info.loc[parity_info['Parity'] == 'Parous', 'Sample'].unique().tolist()
g_samples = parity_info.loc[parity_info['Parity'] == 'G', 'Sample'].unique().tolist()


# =========================
# 6. Final harmonisation function
# =========================
def process_sample(adata):
    adata = adata.copy()

    # Age numeric
    adata.obs['Age'] = pd.to_numeric(adata.obs['Age'], errors='coerce')

    # Add age for tumour samples if needed
    adata.obs.loc[
        adata.obs['Sample'].isin(['WKAL8.2j_tumour_adjacent', 'WKAL8.2j_tumour']),
        'Age'
    ] = 47.0

    # Parity
    adata.obs['Parity'] = 'NP'
    adata.obs.loc[adata.obs['Sample'].isin(parous_samples), 'Parity'] = 'Parous'
    adata.obs.loc[adata.obs['Sample'].isin(g_samples), 'Parity'] = 'G'

    # Tumour
    adata.obs['Tumour'] = 'No'
    adata.obs.loc[adata.obs['Sample'].isin(tumour_samples), 'Tumour'] = 'Tumour'
    adata.obs.loc[adata.obs['Sample'].eq('WKAL8.2j_tumour'), 'Tumour'] = 'Tumour'

    # Genotype
    adata.obs['Genotype'] = adata.obs['Genotype'].astype(str).replace({'WKBR': 'BRCA1'})

    # Age group
    adata.obs['Age_group'] = pd.Series(pd.NA, index=adata.obs.index, dtype='object')
    adata.obs.loc[adata.obs['Age'] < 70, 'Age_group'] = 'Young'
    adata.obs.loc[adata.obs['Age'] >= 70, 'Age_group'] = 'Old'

    # Group
    adata.obs['Group'] = ''

    is_seno = adata.obs['Dataset'].eq('Senolytics')
    is_aging = adata.obs['Dataset'].eq('Aging_atlas')
    is_brca = adata.obs['Dataset'].eq('BRCA2_PALB2')
    is_wkaa = adata.obs['Dataset'].eq('WKAA_First')

    adata.obs.loc[is_seno, 'Group'] = adata.obs.loc[is_seno, 'Treatment'].astype(str)

    adata.obs.loc[is_aging, 'Group'] = (
        adata.obs.loc[is_aging, 'Age_group'].astype(str) + '_' +
        adata.obs.loc[is_aging, 'Parity'].astype(str) + '_' +
        adata.obs.loc[is_aging, 'Genotype'].astype(str)
    )

    adata.obs.loc[is_brca, 'Group'] = adata.obs.loc[is_brca, 'Genotype'].astype(str)
    adata.obs.loc[is_wkaa, 'Group'] = 'WKAA_First'

    group_order = [
        'WT',
        'BRCA2',
        'PALB2',
        'young_NP_BRCA1',
        'Seno_ABT737_control',
        'Seno_ABT737',
        'WKAA_First',
        'young_NP_WT',
        'old_NP_WT',
        'old_Parous_WT',
        'young_G_WT'
    ]

    adata = safe_order_group_categories(adata, group_order)
    return adata


Seno = process_sample(Seno)
WKAA_adata = process_sample(WKAA_adata)
aging_adata = process_sample(aging_adata)
brca_palb2 = process_sample(brca_palb2)

Seno.write_h5ad(f'{OUTDIR}/Seno_adata.h5ad')
WKAA_adata.write_h5ad(f'{OUTDIR}/WKAA_adata.h5ad')
aging_adata.write_h5ad(f'{OUTDIR}/aging_adata.h5ad')
brca_palb2.write_h5ad(f'{OUTDIR}/brca_palb2_adata.h5ad')