import os
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# =============================================================================
# Configuration
# =============================================================================
adata_file = "../integration/QC_files/data.h5ad" 
output_dir = "../integration/QC_results/"
os.makedirs(output_dir, exist_ok=True)

# QC thresholds
MIN_GENES = 500
MIN_COUNTS = 1000
MAX_PCT_MT = 10
MIN_CELLS_PER_GENE = 3

# =============================================================================
# Load data
# =============================================================================
adata = sc.read_h5ad(adata_file)
print(f"Data loaded from {adata_file}")
print(f"  Shape: {adata.shape[0]} cells x {adata.shape[1]} genes")

# Ensure unique barcodes and gene names
adata.obs_names_make_unique()
adata.var_names_make_unique()

# =============================================================================
# Preserve raw counts before any processing
# =============================================================================
adata.layers["counts"] = adata.X.copy()

# =============================================================================
# Annotate gene types
# =============================================================================
adata.var["symbol"] = adata.var.index

# Mitochondrial genes
adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")

# Ribosomal genes
adata.var["ribo"] = adata.var_names.str.upper().str.startswith(("RPS", "RPL"))

# Hemoglobin genes (HBA, HBB, etc., but not HBP)
adata.var["hb"] = adata.var_names.str.upper().str.contains("^HB[AB]", regex=True)

print(f"  Mitochondrial genes: {adata.var['mt'].sum()}")
print(f"  Ribosomal genes: {adata.var['ribo'].sum()}")
print(f"  Hemoglobin genes: {adata.var['hb'].sum()}")

# =============================================================================
# Compute QC metrics
# =============================================================================
print("\nCalculating QC metrics...")
sc.pp.calculate_qc_metrics(
    adata, 
    qc_vars=["mt", "ribo", "hb"], 
    percent_top=[20], 
    log1p=False, 
    inplace=True
)

# Rename columns to standard names
if "total_counts" in adata.obs.columns and "n_counts" not in adata.obs.columns:
    adata.obs.rename(columns={"total_counts": "n_counts"}, inplace=True)

if "n_genes_by_counts" in adata.obs.columns and "n_genes" not in adata.obs.columns:
    adata.obs.rename(columns={"n_genes_by_counts": "n_genes"}, inplace=True)

# Compute counts per gene ratio
adata.obs["counts_per_gene"] = adata.obs["n_counts"] / adata.obs["n_genes"]

# Verify QC metrics
print("\nQC metrics summary:")
print(adata.obs[["n_counts", "n_genes", "pct_counts_mt"]].describe())

# =============================================================================
# Verify required columns exist
# =============================================================================
# Check if mouse_id exists for per-sample plots
has_mouse_id = "mouse_id" in adata.obs.columns
if not has_mouse_id:
    print("\nWARNING: 'mouse_id' column not found. Per-sample plots will be skipped.")
    # Try to find alternative sample column
    for col in ["Sample", "sample", "batch", "Batch"]:
        if col in adata.obs.columns:
            adata.obs["mouse_id"] = adata.obs[col]
            has_mouse_id = True
            print(f"  Using '{col}' as sample identifier.")
            break

# =============================================================================
# Plot QC metrics BEFORE filtering
# =============================================================================
print("\nGenerating pre-QC plots...")

# Cell counts per sample
if has_mouse_id:
    plt.figure(figsize=(12, 6))
    sample_order = adata.obs["mouse_id"].value_counts().index.tolist()
    sns.boxplot(
        x="mouse_id", 
        y=np.log1p(adata.obs["n_counts"]), 
        data=adata.obs, 
        order=sample_order,
        showfliers=False
    )
    plt.xticks(rotation=90)
    plt.xlabel("Sample")
    plt.ylabel("Log Total Counts")
    plt.title("Log Cell Counts Per Sample (Before QC)")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "cell_counts_per_sample_before_qc.png"), dpi=150)
    plt.close()

# Total counts histogram
plt.figure(figsize=(8, 6))
sns.histplot(adata.obs["n_counts"], bins=100, kde=True)
plt.xlabel("Total Counts")
plt.ylabel("Frequency")
plt.title("Total Counts Distribution (Before QC)")
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "total_counts_histogram_before.png"), dpi=150)
plt.close()

# Violin plots for QC metrics
sc.pl.violin(
    adata, 
    keys=["pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"], 
    jitter=0.4, 
    multi_panel=True, 
    show=False
)
plt.savefig(os.path.join(output_dir, "qc_violin_before.png"), bbox_inches="tight", dpi=150)
plt.close()

# Scatter plot: counts vs genes colored by MT%
plt.figure(figsize=(8, 6))
scatter = plt.scatter(
    np.log1p(adata.obs["n_counts"]), 
    np.log1p(adata.obs["n_genes"]), 
    c=adata.obs["pct_counts_mt"], 
    cmap="viridis", 
    alpha=0.5,
    s=1
)
plt.colorbar(scatter, label="MT%")
plt.xlabel("Log Total Counts")
plt.ylabel("Log Genes Detected")
plt.title("Counts vs Genes (Before QC)")
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "qc_scatter_log_before.png"), dpi=150)
plt.close()

# Counts per gene distribution
plt.figure(figsize=(8, 6))
plt.hist(adata.obs["counts_per_gene"], bins=100, color="skyblue", edgecolor="black", linewidth=0.5)
plt.xlabel("Counts per Gene")
plt.ylabel("Frequency")
plt.title("Counts per Gene Distribution (Before QC)")
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "counts_per_gene_before.png"), dpi=150)
plt.close()

# Mitochondrial percentage distribution
plt.figure(figsize=(8, 6))
plt.hist(adata.obs["pct_counts_mt"], bins=100, color="salmon", edgecolor="black", linewidth=0.5)
plt.axvline(x=MAX_PCT_MT, color="red", linestyle="--", label=f"Threshold: {MAX_PCT_MT}%")
plt.xlabel("Mitochondrial Percentage")
plt.ylabel("Frequency")
plt.title("Mitochondrial Percentage Distribution")
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "mitochondrial_percentage_distribution.png"), dpi=150)
plt.close()

# =============================================================================
# Define QC pass/fail
# =============================================================================
qc_pass = (
    (adata.obs["n_genes"] >= MIN_GENES) & 
    (adata.obs["pct_counts_mt"] <= MAX_PCT_MT)
)

# Store QC status (True = passed, False = failed)
adata.obs["QC_passed"] = qc_pass.astype(int)
adata.obs["QC_failed"] = (~qc_pass).astype(int)

# =============================================================================
# Apply filtering
# =============================================================================
print("\n" + "=" * 50)
print("FILTERING")
print("=" * 50)

before_cells = adata.n_obs
before_genes = adata.n_vars

# Filter cells
adata_filtered = adata[qc_pass, :].copy()

# Filter genes (remove genes expressed in very few cells)
sc.pp.filter_genes(adata_filtered, min_cells=MIN_CELLS_PER_GENE)

after_cells = adata_filtered.n_obs
after_genes = adata_filtered.n_vars

print(f"Cells: {before_cells} -> {after_cells} ({before_cells - after_cells} removed, {100 * (before_cells - after_cells) / before_cells:.1f}%)")
print(f"Genes: {before_genes} -> {after_genes} ({before_genes - after_genes} removed)")

# Breakdown of filtering reasons
n_low_genes = (adata.obs["n_genes"] < MIN_GENES).sum()
n_high_mt = (adata.obs["pct_counts_mt"] > MAX_PCT_MT).sum()

print(f"\nFiltering breakdown:")
print(f"  Low gene count (< {MIN_GENES}): {n_low_genes}")
print(f"  High MT% (> {MAX_PCT_MT}%): {n_high_mt}")

# =============================================================================
# Plot QC metrics AFTER filtering
# =============================================================================
print("\nGenerating post-QC plots...")

# Recalculate counts_per_gene for filtered data
adata_filtered.obs["counts_per_gene"] = adata_filtered.obs["n_counts"] / adata_filtered.obs["n_genes"]

# Cell counts per sample 
if has_mouse_id:
    plt.figure(figsize=(12, 6))
    sample_order = adata_filtered.obs["mouse_id"].value_counts().index.tolist()
    sns.boxplot(
        x="mouse_id", 
        y=np.log1p(adata_filtered.obs["n_counts"]), 
        data=adata_filtered.obs,
        order=sample_order,
        showfliers=False
    )
    plt.xticks(rotation=90)
    plt.xlabel("Sample")
    plt.ylabel("Log Total Counts")
    plt.title("Log Cell Counts Per Sample (After QC)")
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "cell_counts_per_sample_after_qc.png"), dpi=150)
    plt.close()

# Total counts histogram
plt.figure(figsize=(8, 6))
sns.histplot(adata_filtered.obs["n_counts"], bins=100, kde=True)
plt.xlabel("Total Counts")
plt.ylabel("Frequency")
plt.title("Total Counts Distribution (After QC)")
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "total_counts_histogram_after.png"), dpi=150)
plt.close()

# Violin plots
sc.pl.violin(
    adata_filtered, 
    keys=["pct_counts_mt", "pct_counts_ribo", "pct_counts_hb"], 
    jitter=0.4, 
    multi_panel=True, 
    show=False
)
plt.savefig(os.path.join(output_dir, "qc_violin_after.png"), bbox_inches="tight", dpi=150)
plt.close()

# Scatter plot
plt.figure(figsize=(8, 6))
scatter = plt.scatter(
    np.log1p(adata_filtered.obs["n_counts"]), 
    np.log1p(adata_filtered.obs["n_genes"]), 
    c=adata_filtered.obs["pct_counts_mt"], 
    cmap="viridis", 
    alpha=0.5,
    s=1
)
plt.colorbar(scatter, label="MT%")
plt.xlabel("Log Total Counts")
plt.ylabel("Log Genes Detected")
plt.title("Counts vs Genes (After QC)")
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "qc_scatter_log_after.png"), dpi=150)
plt.close()

# Counts per gene distribution
plt.figure(figsize=(8, 6))
plt.hist(adata_filtered.obs["counts_per_gene"], bins=100, color="skyblue", edgecolor="black", linewidth=0.5)
plt.xlabel("Counts per Gene")
plt.ylabel("Frequency")
plt.title("Counts per Gene Distribution (After QC)")
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "counts_per_gene_after.png"), dpi=150)
plt.close()

# =============================================================================
# UMAP visualization of QC status (on full data)
# =============================================================================
print("\nComputing UMAP for QC visualization...")

# Work on a copy for UMAP to preserve original
adata_umap = adata.copy()

# Normalize and log-transform (using the copy)
sc.pp.normalize_total(adata_umap, target_sum=1e4)
sc.pp.log1p(adata_umap)

# Find highly variable genes
sc.pp.highly_variable_genes(adata_umap, n_top_genes=2000)

# PCA, neighbors, UMAP
sc.pp.pca(adata_umap, n_comps=50, random_state=42)
sc.pp.neighbors(adata_umap, n_neighbors=15, n_pcs=50, random_state=42)
sc.tl.umap(adata_umap, random_state=42)

# Plot UMAPs
sc.settings.figdir = output_dir

sc.pl.umap(
    adata_umap, 
    color=["QC_passed"], 
    title="QC Status (1=Passed, 0=Failed)",
    save="_QC_status.png", 
    show=False
)

sc.pl.umap(
    adata_umap, 
    color=["n_counts", "n_genes", "pct_counts_mt", "pct_counts_ribo", "pct_counts_hb", "counts_per_gene"],
    save="_QC_metrics.png", 
    show=False
)

# =============================================================================
# Save outputs
# =============================================================================
print("\nSaving outputs...")

# Save filtered data
filtered_path = os.path.join(output_dir, "adata_qc_filtered.h5ad")
adata_filtered.write(filtered_path)
print(f"  Filtered data: {filtered_path}")

# Save QC summary
qc_summary = pd.DataFrame({
    "metric": [
        "total_cells_before",
        "total_cells_after", 
        "cells_removed",
        "pct_cells_removed",
        "total_genes_before",
        "total_genes_after",
        "genes_removed",
        "mean_counts_after",
        "median_counts_after",
        "mean_genes_after",
        "median_genes_after",
        "mean_pct_mt_after",
        "median_pct_mt_after",
        "mean_pct_ribo_after",
        "median_pct_ribo_after",
    ],
    "value": [
        before_cells,
        after_cells,
        before_cells - after_cells,
        100 * (before_cells - after_cells) / before_cells,
        before_genes,
        after_genes,
        before_genes - after_genes,
        adata_filtered.obs["n_counts"].mean(),
        adata_filtered.obs["n_counts"].median(),
        adata_filtered.obs["n_genes"].mean(),
        adata_filtered.obs["n_genes"].median(),
        adata_filtered.obs["pct_counts_mt"].mean(),
        adata_filtered.obs["pct_counts_mt"].median(),
        adata_filtered.obs["pct_counts_ribo"].mean(),
        adata_filtered.obs["pct_counts_ribo"].median(),
    ]
})

summary_path = os.path.join(output_dir, "qc_summary.csv")
qc_summary.to_csv(summary_path, index=False)
print(f"  QC summary: {summary_path}")

# Per-sample summary
if has_mouse_id:
    per_sample = adata_filtered.obs.groupby("mouse_id").agg(
        n_cells=("n_counts", "count"),
        mean_counts=("n_counts", "mean"),
        mean_genes=("n_genes", "mean"),
        mean_pct_mt=("pct_counts_mt", "mean"),
    ).reset_index()
    
    per_sample_path = os.path.join(output_dir, "qc_per_sample.csv")
    per_sample.to_csv(per_sample_path, index=False)
    print(f"  Per-sample summary: {per_sample_path}")

print("\nQC complete! Data ready for doublet detection.")
