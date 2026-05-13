import os
import scanpy as sc
import numpy as np
import pandas as pd
import scrublet as scr
import matplotlib.pyplot as plt
import seaborn as sns

# =============================================================================
# Configuration
# =============================================================================
input_file = '../adata/data.h5ad'  
output_dir = "../integration/adatas"
doublet_dir = os.path.join(output_dir, "doublet_results")
os.makedirs(doublet_dir, exist_ok=True)

# Scrublet parameters
DEFAULT_THRESHOLD = 0.33
RANDOM_STATE = 42

# =============================================================================
# Load data
# =============================================================================
adata = sc.read_h5ad(input_file)
print(f"Data loaded from {input_file}")
print(f"  Shape: {adata.shape[0]} cells x {adata.shape[1]} genes")

# =============================================================================
# Verify raw counts (Scrublet requires raw counts, not normalized data)
# =============================================================================
# Check if data looks like raw counts
sample_values = adata.X[:100, :100].toarray() if hasattr(adata.X, 'toarray') else adata.X[:100, :100]
if np.any(sample_values < 0) or np.any(sample_values % 1 != 0):
    print("WARNING: Data may not be raw counts. Scrublet expects integer count data.")
    print("  If you have a 'counts' layer, consider using: adata.X = adata.layers['counts']")

# =============================================================================
# Compute QC metrics if missing
# =============================================================================
required_qc_cols = ["n_counts", "n_genes", "pct_counts_mt"]
missing_cols = [col for col in required_qc_cols if col not in adata.obs.columns]

if missing_cols:
    print(f"Missing QC metrics: {missing_cols}. Recomputing...")
    
    # Ensure mitochondrial genes are correctly annotated
    if "mt" not in adata.var.columns:
        adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    
    # Compute QC metrics
    sc.pp.calculate_qc_metrics(
        adata, 
        qc_vars=["mt"], 
        percent_top=None, 
        log1p=False, 
        inplace=True
    )
    
    # Rename columns to match expected format
    if "total_counts" in adata.obs.columns and "n_counts" not in adata.obs.columns:
        adata.obs.rename(columns={"total_counts": "n_counts"}, inplace=True)
    if "n_genes_by_counts" in adata.obs.columns and "n_genes" not in adata.obs.columns:
        adata.obs.rename(columns={"n_genes_by_counts": "n_genes"}, inplace=True)
    
    print("QC metrics successfully recomputed.")

# =============================================================================
# Initialize doublet detection columns
# =============================================================================
adata.obs["doublet_score"] = np.zeros(adata.n_obs, dtype=np.float32)
adata.obs["predicted_doublet"] = np.full(adata.n_obs, False, dtype=bool)
adata.obs["doublet_threshold"] = np.full(adata.n_obs, np.nan, dtype=np.float32)

# =============================================================================
# Run Scrublet per sample
# =============================================================================
samples = adata.obs["Sample"].unique()
n_samples = len(samples)

print(f"\nRunning Scrublet on {n_samples} samples...")

for i, sample in enumerate(samples, 1):
    print(f"[{i}/{n_samples}] Processing sample: {sample}")
    
    # Get mask and subset data
    mask = adata.obs["Sample"] == sample
    n_cells = mask.sum()
    print(f"  Cells: {n_cells}")
    
    if n_cells < 50:
        print(f"  Skipping: too few cells for reliable doublet detection")
        continue
    
    # Create a copy to avoid view issues
    lane_adata = adata[mask, :].copy()
    
    # Run Scrublet with explicit random state
    try:
        scrub = scr.Scrublet(lane_adata.X, random_state=RANDOM_STATE)
        doublet_scores, predicted_doublets = scrub.scrub_doublets()
    except Exception as e:
        print(f"  ERROR running Scrublet: {e}")
        continue
    
    # Handle threshold
    if scrub.threshold_ is None:
        print(f"  No automatic threshold found, using default: {DEFAULT_THRESHOLD}")
        predicted_doublets = scrub.call_doublets(threshold=DEFAULT_THRESHOLD)
        threshold_value = DEFAULT_THRESHOLD
    else:
        threshold_value = scrub.threshold_
        print(f"  Threshold: {threshold_value:.3f}")
    
    # Get indices for assignment
    indices = adata.obs.index[mask]
    
    # Store results
    adata.obs.loc[indices, "doublet_score"] = doublet_scores.astype(np.float32)
    adata.obs.loc[indices, "predicted_doublet"] = predicted_doublets
    adata.obs.loc[indices, "doublet_threshold"] = threshold_value
    
    # Report doublet count
    n_doublets = predicted_doublets.sum()
    pct_doublets = 100 * n_doublets / n_cells
    print(f"  Doublets: {n_doublets} ({pct_doublets:.1f}%)")

# =============================================================================
# Create doublet status column
# =============================================================================
adata.obs["doublet_status"] = np.where(
    adata.obs["predicted_doublet"], 
    "Doublet", 
    "Singlet"
)
adata.obs["doublet_status"] = adata.obs["doublet_status"].astype("category")

# =============================================================================
# Summary statistics
# =============================================================================
print("\n" + "=" * 50)
print("SUMMARY")
print("=" * 50)
total_doublets = adata.obs["predicted_doublet"].sum()
total_cells = adata.n_obs
print(f"Total cells: {total_cells}")
print(f"Total doublets: {total_doublets} ({100 * total_doublets / total_cells:.1f}%)")
print(f"Total singlets: {total_cells - total_doublets}")

# Per-sample summary
summary_df = adata.obs.groupby("Sample").agg(
    n_cells=("predicted_doublet", "count"),
    n_doublets=("predicted_doublet", "sum"),
    mean_score=("doublet_score", "mean"),
    threshold=("doublet_threshold", "first")
).reset_index()
summary_df["pct_doublets"] = 100 * summary_df["n_doublets"] / summary_df["n_cells"]
print("\nPer-sample summary:")
print(summary_df.to_string(index=False))

# =============================================================================
# Save doublet scores
# =============================================================================
doublet_scores_path = os.path.join(doublet_dir, "doublet_scores.csv")
adata.obs[["Sample", "doublet_score", "doublet_status", "doublet_threshold"]].to_csv(
    doublet_scores_path
)
print(f"\nDoublet scores saved to {doublet_scores_path}")

# Save summary
summary_path = os.path.join(doublet_dir, "doublet_summary.csv")
summary_df.to_csv(summary_path, index=False)
print(f"Summary saved to {summary_path}")

# =============================================================================
# Visualizations
# =============================================================================
print("\nGenerating visualizations...")

# 1. Per-sample histograms
for sample in samples:
    lane_data = adata.obs.loc[adata.obs["Sample"] == sample, ["doublet_score", "doublet_threshold"]]
    
    if lane_data["doublet_score"].isna().all():
        continue
    
    plt.figure(figsize=(8, 6))
    sns.histplot(lane_data["doublet_score"].dropna(), bins=50, kde=True)
    
    # Plot threshold line
    threshold = lane_data["doublet_threshold"].dropna().unique()
    if len(threshold) == 1 and not np.isnan(threshold[0]):
        plt.axvline(
            x=threshold[0], 
            color="red", 
            linestyle="--", 
            label=f"Threshold: {threshold[0]:.2f}"
        )
    
    plt.xlabel("Doublet Score")
    plt.ylabel("Frequency")
    plt.title(f"Scrublet Doublet Score Distribution ({sample})")
    plt.legend()
    
    # Clean filename
    safe_sample_name = str(sample).replace("/", "_").replace("\\", "_")
    plt.savefig(
        os.path.join(doublet_dir, f"doublet_score_distribution_{safe_sample_name}.png"),
        bbox_inches="tight",
        dpi=150
    )
    plt.close()

# 2. Combined overview plot
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Histogram of all doublet scores
ax1 = axes[0]
sns.histplot(adata.obs["doublet_score"].dropna(), bins=50, kde=True, ax=ax1)
ax1.set_xlabel("Doublet Score")
ax1.set_ylabel("Frequency")
ax1.set_title("Overall Doublet Score Distribution")

# Doublet percentage per sample
ax2 = axes[1]
summary_df_sorted = summary_df.sort_values("pct_doublets", ascending=True)
ax2.barh(summary_df_sorted["Sample"].astype(str), summary_df_sorted["pct_doublets"])
ax2.set_xlabel("Doublet Percentage (%)")
ax2.set_ylabel("Sample")
ax2.set_title("Doublet Rate per Sample")

# Doublet score vs n_counts (if available)
ax3 = axes[2]
if "n_counts" in adata.obs.columns:
    # Subsample for plotting if too many cells
    if adata.n_obs > 10000:
        plot_idx = np.random.choice(adata.n_obs, 10000, replace=False)
        plot_data = adata.obs.iloc[plot_idx]
    else:
        plot_data = adata.obs
    
    scatter = ax3.scatter(
        plot_data["n_counts"],
        plot_data["doublet_score"],
        c=plot_data["predicted_doublet"].astype(int),
        cmap="coolwarm",
        alpha=0.3,
        s=1
    )
    ax3.set_xlabel("Total Counts")
    ax3.set_ylabel("Doublet Score")
    ax3.set_title("Doublet Score vs Total Counts")
else:
    ax3.text(0.5, 0.5, "n_counts not available", ha="center", va="center", transform=ax3.transAxes)
    ax3.set_title("Doublet Score vs Total Counts")

plt.tight_layout()
plt.savefig(os.path.join(doublet_dir, "doublet_overview.png"), bbox_inches="tight", dpi=150)
plt.close()

print(f"Visualizations saved to {doublet_dir}")

# =============================================================================
# Save annotated adata
# =============================================================================
updated_adata_path = os.path.join(doublet_dir, "adata_with_doublet_scores.h5ad")
adata.write(updated_adata_path)
print(f"\nAnnotated AnnData saved to {updated_adata_path}")

print("\nDone!")
