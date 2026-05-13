import os
import random

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import scvi
import torch


# ----------------------------
# Settings
# ----------------------------
SEED = 42
np.random.seed(SEED)
random.seed(SEED)
torch.manual_seed(SEED)
if torch.cuda.is_available():
    torch.cuda.manual_seed_all(SEED)

scvi.settings.seed = SEED
sc.set_figure_params(figsize=(6, 6))

accelerator = "gpu" if torch.cuda.is_available() else "cpu"

BASE_DIR = "../scvi_int/all/output"
DATA_DIR = "../data"
OLD_META_PATH = "../data/old_meta/meta.csv"

INPUT_H5AD = os.path.join(DATA_DIR, "adata_full_no_varm.h5ad")
PREP_H5AD = os.path.join(DATA_DIR, "adata_full_pre_scvi.h5ad")

SAVE_DIR = BASE_DIR
LOSS_DIR = os.path.join(SAVE_DIR, "loss")
UMAP_DIR = os.path.join(SAVE_DIR, "figures")
ADATA_DIR = os.path.join(SAVE_DIR, "adata")
MODEL_DIR = os.path.join(SAVE_DIR, "model")

FINAL_ADATA_PATH = os.path.join(ADATA_DIR, "adata_full_scvi.h5ad")
FINAL_MODEL_PATH = os.path.join(MODEL_DIR, "adata_full_scvi_model")
MOUSE_INFO_PATH = os.path.join(SAVE_DIR, "before_scvi_mouse_info.csv")
GENE_MAP_PATH = os.path.join(SAVE_DIR, "ID_to_Symbol_genes.csv")

for p in [LOSS_DIR, UMAP_DIR, ADATA_DIR, MODEL_DIR]:
    os.makedirs(p, exist_ok=True)

sc.settings.figdir = UMAP_DIR


# ----------------------------
# Helper functions
# ----------------------------
def run_scvi(
    adata,
    batch_key="Batch",
    n_dims=20,
    n_layers=2,
    max_epochs=600,
    resume_path=None,
    save_path=None,
):
    if "highly_variable" not in adata.var.columns:
        raise ValueError("'highly_variable' not found in adata.var. Run HVG selection first.")

    if "raw" not in adata.layers:
        raise ValueError("'raw' layer not found. scVI should be trained on raw counts.")

    train_adata = adata[:, adata.var["highly_variable"]].copy()
    train_adata.X = train_adata.layers["raw"].copy()

    scvi.model.SCVI.setup_anndata(train_adata, batch_key=batch_key)

    model_params = {
        "n_latent": n_dims,
        "gene_likelihood": "nb",
        "use_layer_norm": "both",
        "use_batch_norm": "none",
        "encode_covariates": True,
        "dropout_rate": 0.2,
        "n_layers": n_layers,
    }

    if resume_path is not None:
        print(f"Resuming training from {resume_path}")
        vae = scvi.model.SCVI.load(resume_path, adata=train_adata)
    else:
        vae = scvi.model.SCVI(train_adata, **model_params)

    vae.train(
        early_stopping=True,
        train_size=0.9,
        early_stopping_patience=45,
        max_epochs=max_epochs,
        batch_size=1024,
        accelerator=accelerator,
    )

    latent = vae.get_latent_representation(train_adata)
    if latent.shape[0] != adata.n_obs:
        raise ValueError("Latent representation row count does not match number of cells.")

    adata.obsm["X_scVI"] = latent

    if save_path is not None:
        vae.save(save_path, overwrite=True)

    return adata, vae


def save_training_curves(vae, outdir):
    plt.clf()
    plt.plot(vae.history["elbo_train"], label="train", color="blue")
    plt.plot(vae.history["elbo_validation"], label="validation", color="red")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "elbo.png"), dpi=300)

    plt.clf()
    plt.plot(vae.history["reconstruction_loss_train"], label="train", color="blue")
    plt.plot(vae.history["reconstruction_loss_validation"], label="validation", color="red")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "reconstruction.png"), dpi=300)


def clear_slots(adata):
    adata.obsm.clear()
    adata.obsp.clear()
    adata.varm.clear()
    return adata


# ----------------------------
# Load data
# ----------------------------
adata_full = ad.read_h5ad(INPUT_H5AD)
adata_full = adata_full.copy()

adata_full.obs_names_make_unique()
adata_full.var_names_make_unique()

if "raw" not in adata_full.layers:
    raise ValueError("Expected raw counts in adata_full.layers['raw'].")


# ----------------------------
# Fix metadata
# ----------------------------
adata_full.obs["Genotype"] = adata_full.obs["Genotype"].astype(str)
adata_full.obs.loc[
    adata_full.obs["Mouse_ID"].isin(["WKAL8.2j_tumour"]),
    "Genotype"
] = "BRCA2_Tumour"
adata_full.obs.loc[
    adata_full.obs["Mouse_ID"].isin(["WKBR75.4bTM", "WKBR61.4dTM"]),
    "Genotype"
] = "BRCA1_Tumour"
adata_full.obs.loc[
    adata_full.obs["Mouse_ID"].isin(["WKAL8.2j_tumour_adjacent", "WKAL8.2j_tumour"]),
    "Age"
] = 47.0

adata_full.obs["Parity"] = "NP"
adata_full.obs.loc[
    adata_full.obs["Mouse_ID"].isin(["SIGAD10", "SIGAB10", "SIGAC10", "SIGAA10"]),
    "Parity"
] = "Parous"
adata_full.obs.loc[
    adata_full.obs["Mouse_ID"].isin(["4.5dG_3", "4.5dG_2", "4.5dG_1"]),
    "Parity"
] = "4.5dG"
adata_full.obs.loc[
    adata_full.obs["Mouse_ID"].isin(["9.5dG_3", "9.5dG_1", "9.5dG_2"]),
    "Parity"
] = "9.5dG"
adata_full.obs.loc[
    adata_full.obs["Mouse_ID"].isin(["14.5dG_3", "14.5dG_1", "14.5dG_2"]),
    "Parity"
] = "14.5dG"

mouse_info = adata_full.obs[
    ["Mouse_ID", "Batch", "Dataset", "Age", "Genotype", "Treatment", "Parity"]
].drop_duplicates()
mouse_info.to_csv(MOUSE_INFO_PATH, index=False)


# ----------------------------
# Merge old metadata
# ----------------------------
# old_meta = pd.read_csv(OLD_META_PATH, index_col=0)
# old_meta = old_meta[["Condition", "subsets_Mito_percent"]]
# adata_full.obs = adata_full.obs.join(old_meta, how="left")


# ----------------------------
# QC filtering
# ----------------------------
if "pct_counts_mt" not in adata_full.obs.columns:
    raise ValueError("'pct_counts_mt' not found in adata_full.obs.")

cells_to_remove = adata_full.obs.index[adata_full.obs["pct_counts_mt"] > 10]
print(f"Number of cells to remove: {len(cells_to_remove)}")
adata_full = adata_full[~adata_full.obs.index.isin(cells_to_remove)].copy()

# ----------------------------
# Normalization + HVGs for plotting / downstream
# ----------------------------
adata_full.X = adata_full.layers["raw"].copy()
sc.pp.normalize_total(adata_full, target_sum=1e4)
sc.pp.log1p(adata_full)
sc.pp.highly_variable_genes(adata_full, n_top_genes=5000)
adata_full.layers["log1p_norm"] = adata_full.X.copy()

adata_full.write_h5ad(PREP_H5AD)


# ----------------------------
# UMAP before scVI
# ----------------------------
sc.pp.pca(adata_full, n_comps=50, use_highly_variable=True, svd_solver="arpack")
sc.pp.neighbors(adata_full, n_neighbors=15, n_pcs=50)
sc.tl.umap(adata_full)
sc.pl.umap(adata_full, color="Batch", save="_before_scvi_batch.png", show=False)


# ----------------------------
# Run scVI
# ----------------------------
adata_full = ad.read_h5ad(PREP_H5AD)
adata_full, vae = run_scvi(
    adata_full,
    batch_key="Batch",
    n_dims=20,
    n_layers=2,
    max_epochs=600,
    save_path=FINAL_MODEL_PATH,
)

save_training_curves(vae, LOSS_DIR)


# ----------------------------
# UMAP after scVI
# ----------------------------
sc.pp.neighbors(adata_full, use_rep="X_scVI", n_neighbors=15)
sc.tl.umap(adata_full)
sc.pl.umap(adata_full, color="Batch", save="_adata_full_batch.png", show=False)


# ----------------------------
# Save main outputs
# ----------------------------
adata_full.write_h5ad(FINAL_ADATA_PATH)


# ----------------------------
# Reload and add gene symbols
# ----------------------------
adata_full = ad.read_h5ad(FINAL_ADATA_PATH)

if os.path.exists(GENE_MAP_PATH):
    map_df = pd.read_csv(GENE_MAP_PATH, index_col=0)
    adata_full.var = adata_full.var.merge(map_df, left_index=True, right_index=True, how="left")

    if "Symbol" in adata_full.var.columns:
        adata_full.var["Symbol"] = adata_full.var["Symbol"].fillna(adata_full.var_names)
        adata_full.var_names = adata_full.var["Symbol"].astype(str)
        adata_full.var_names_make_unique()


# ----------------------------
# Leiden clustering
# ----------------------------
leiden_resolutions = [0.05, 0.1, 0.2, 0.3, 0.4, 0.5]
for res in leiden_resolutions:
    key = "leiden_" + str(res).replace(".", "_")
    sc.tl.leiden(adata_full, resolution=res, key_added=key)
    sc.pl.umap(adata_full, color=key, save=f"_{key}.png", show=False)


# ----------------------------
# Marker dotplots
# ----------------------------
markers_dict = {
    "Epithelial": ["Krt8", "Sfn", "Krt18", "Epcam"],
    "Fibroblast": ["Col3a1", "Dcn", "Col4a1"],
    "Immune": ["Ptprc", "Cd52"],
    "Myeloid": ["Tyrobp", "H2-Aa"],
    "Lymphoid": ["Hcst", "Cd2", "Cd3d"],
    "Stroma": ["Pecam1", "Cdh5", "Pdgfra", "Rgs5"],
}

if "log1p_norm" not in adata_full.layers:
    raise ValueError("'log1p_norm' layer not found.")

adata_full.X = adata_full.layers["log1p_norm"].copy()

for res in leiden_resolutions:
    key = "leiden_" + str(res).replace(".", "_")
    print(f"Running dendrogram and dotplot for resolution {res}")
    sc.tl.dendrogram(adata_full, groupby=key, use_rep="X_scVI")
    sc.pl.dotplot(
        adata_full,
        groupby=key,
        var_names=markers_dict,
        dendrogram=True,
        save=f"_dotplot_{key}.png",
        show=False,
    )

# ----------------------------
# Process obs adding final group column
# ----------------------------

match_dict = {
        'old_Parous_WT_Aging_atlas_None': 'Old_WT_Parous',
        'old_NP_WT_Aging_atlas_None': 'Old_WT_NP',
        'young_NP_WT_Aging_atlas_None': 'Young_WT_NP',
        'young_NP_BRCA1_Tumour_Aging_atlas_None': 'BRCA1_Tumour',
        'young_NP_BRCA1_Aging_atlas_None': 'BRCA1',
        'young_9.5dG_WT_Aging_atlas_None': 'Young_WT_9.5dG',
        'young_4.5dG_WT_Aging_atlas_None': 'Young_WT_4.5dG',
        'young_14.5dG_WT_Aging_atlas_None': 'Young_WT_14.5dG',
        'young_NP_PALB2_BRCA2_PALB2_None': 'PALB2',
        'young_NP_BRCA2_BRCA2_PALB2_None': 'BRCA2',
        'young_NP_BRCA2_Tumour_BRCA2_PALB2_None': 'BRCA2_Tumour',
        'young_NP_WT_BRCA2_PALB2_None': 'Young_WT_NP',
        'young_NP_BRCA1_Senolytics_Seno_ABT737_control': 'ABT737_control',
        'young_NP_BRCA1_Senolytics_Seno_ABT737': 'Treated',
        'young_NP_BRCA1_WKAA_First_None': 'BRCA1'
    }

def process_adata(adata):
    adata.obs['age_group'] = np.where(adata.obs['Age'] < 70, 'young', 'old')
    adata.obs['Genotype_old'] = adata.obs['Genotype'].copy()
    adata.obs['Genotype'] = adata.obs['Genotype'].astype(str)
    adata.obs.loc[adata.obs['Genotype'] == 'WKBR', 'Genotype'] = 'BRCA1'
    adata.obs.loc[adata.obs['Genotype'] == 'WKAA', 'Genotype'] = 'BRCA1'
    adata.obs['Group_long'] = adata.obs['age_group'].astype(str) + '_' + adata.obs['Parity'].astype(str) + '_' + adata.obs['Genotype'].astype(str) + '_' + adata.obs['Dataset'].astype(str) + '_' + adata.obs['Treatment'].astype(str)
    adata.obs['Group'] = adata.obs['Group_long'].map(match_dict)
    return adata

adata_full = process_adata(adata_full)


# ----------------------------
# Manual annotation
# ----------------------------

##The manual annotation was performed first by looking at the expression of known marker genes across clusters, 
##and then by looking at the DE genes for each cluster (using sc.tl.rank_genes_groups with the leiden clusters as groups).

anno_dict = {
    "0": "Stromal",
    "1": "Stromal",
    "2": "Epithelial",
    "3": "Stromal",
    "4": "Lymphoid",
    "5": "Lymphoid",
    "6": "Epithelial",
    "7": "Stromal",
    "8": "Epithelial",
    "9": "Epithelial",
    "10": "Epithelial",
    "11": "Myeloid",
    "12": "Epithelial",
    "13": "Lymphoid",
    "14": "Lymphoid",
    "15": "Stromal",
    "16": "Myeloid",
    "17": "Epithelial",
    "18": "Doublet",
    "19": "Epithelial",
}

adata_full.obs["ct_level1"] = adata_full.obs["leiden_0_2"].astype(str).map(anno_dict)
sc.pl.umap(adata_full, color="ct_level1", save="_umap_ct_level1.png", show=False)

adata_full.write_h5ad(FINAL_ADATA_PATH)


# ----------------------------
# Export major compartments
# ----------------------------
adata_full = ad.read_h5ad(FINAL_ADATA_PATH)

subsets = {
    "lym": "Lymphoid",
    "mye": "Myeloid",
    "epi": "Epithelial",
    "stromal": "Stromal",
}

for folder, label in subsets.items():
    sub = adata_full[adata_full.obs["ct_level1"] == label].copy()
    sub = clear_slots(sub)
    sub.X = sub.layers["raw"].copy()

    subset_dir = os.path.join("../scvi_int", folder)
    os.makedirs(subset_dir, exist_ok=True)
    sub.write_h5ad(os.path.join(subset_dir, f"{folder}_scvi.h5ad"))
