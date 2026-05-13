import os
import re
from pathlib import Path
import numpy as np
import anndata as ad
import scanpy as sc
import scipy.sparse as sp
import matplotlib as mpl
import matplotlib.pyplot as plt


def _clean(name: str) -> str:
    return re.sub(r'[^a-zA-Z0-9._-]+', '_', name).strip('_')

def plot_dotplot(
    adata: ad.AnnData,
    genes: list,
    name: str,
    groupby: str = 'Group',
    ct_key: str = 'ct_level3',
    output_dir: str = '../analysis/senescence/dotplot',
    min_cells_show: int = 1,
    min_cells_reliable: int = 10,
    dot_max: float = 0.5,
    vmax: float = 1.0,
):
    """
    Creates one dotplot per ct_level3, masking expression to zero for groups with
    min_cells_show <= n < min_cells_reliable so they appear as empty dots.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    sc.settings.figdir = str(output_dir) 

    style_rc = {
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial'],
        'axes.edgecolor': 'black',
        'axes.linewidth': 0.8,
        'xtick.major.size': 2,
        'ytick.major.size': 2,
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'figure.dpi': 300,
        'savefig.dpi': 300,
        'axes.titlesize': 10,
        'axes.labelsize': 9,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'legend.fontsize': 8,
    }

    ##load and preprocess data
    adata.X = adata.layers['log1p_norm']
    adata = adata[adata.obs['Group'].isin(['Old_WT_Parous','Old_WT_NP','Young_WT_NP','BRCA1','PALB2','BRCA2'])]
    adata.obs['Group'] = adata.obs['Group'].astype('category').cat.reorder_categories(
        ['Young_WT_NP','Old_WT_NP','Old_WT_Parous','BRCA1','BRCA2','PALB2'])
    adata.obs['ct_level3'] = adata.obs['ct_level3'].astype(str)

    ##check genes and report missing ones
    present_genes = [g for g in genes if g in adata.var_names]
    missing = [g for g in genes if g not in adata.var_names]
    if missing:
        print(f"[{name}] {len(missing)} genes not found and will be skipped: {', '.join(missing[:10])}" +
              (" ..." if len(missing) > 10 else ""))
    if len(present_genes) == 0:
        print(f"[{name}] No genes present in var_names. Skipping.")
        return

    for ct in adata.obs[ct_key].unique().tolist():
        ct_mask = (adata.obs[ct_key] == ct).values
        if ct_mask.sum() == 0:
            continue

        print(f"[{name}] Processing cell type: {ct} (n={ct_mask.sum()})")
        cell_data = adata[ct_mask].copy()  

        if not hasattr(cell_data.obs[groupby].dtype, "categories"):
            cell_data.obs[groupby] = cell_data.obs[groupby].astype('category')
        group_sizes = cell_data.obs[groupby].value_counts().sort_index()
        valid_groups = group_sizes[group_sizes >= min_cells_show].index.tolist()
        if len(valid_groups) <= 1:
            continue

        unreliable_groups = group_sizes[(group_sizes >= min_cells_show) & (group_sizes < min_cells_reliable)].index.tolist()

        if sp.issparse(cell_data.X):
            X = cell_data.X.toarray()
        else:
            X = np.asarray(cell_data.X)

        gene_idx = [cell_data.var_names.get_loc(g) for g in present_genes]

        if len(unreliable_groups) > 0 and len(gene_idx) > 0:
            unreli_mask = cell_data.obs[groupby].isin(unreliable_groups).values
            if unreli_mask.any():
                X[np.where(unreli_mask)[0][:, None], gene_idx] = 0.0

        cell_data.layers["masked"] = sp.csr_matrix(X) if not sp.issparse(cell_data.X) else sp.csr_matrix(X)

        cell_data = cell_data[cell_data.obs[groupby].isin(valid_groups)].copy()

        ct_str = str(ct).replace('/', '_')
        title = ct_str.replace('_', ' ')
        suffix_pdf = f"{_clean(ct_str)}_age_dotplot_{_clean(name)}.pdf"
        suffix_svg = f"{_clean(ct_str)}_age_dotplot_{_clean(name)}.svg"

        ##plot and save
        with mpl.rc_context(style_rc):
            dp = sc.pl.dotplot(
                cell_data,
                var_names=present_genes,
                groupby=groupby,
                use_raw=False,
                layer="masked",
                dot_max=dot_max,
                vmax=vmax,
                standard_scale='var',
                show=False,
                return_fig=True,
                title=title
            )
            dp.savefig(output_dir / suffix_pdf, bbox_inches="tight")
            dp.savefig(output_dir / suffix_svg, bbox_inches="tight")
            plt.close(dp.fig)

    print(f"[{name}] Done. Figures in: {output_dir}")

immune_receptors = ["Pdcd1","Ctla4","Tigit","Lag3","Cd96"]
immune_ligands   = ["Cd274","Pdcd1lg2","Vtcn1","Fgl1","Lgals9"]

epi = ad.read_h5ad('../epi_scvi_post.h5ad')
lym = ad.read_h5ad('../lym_post_scvi.h5ad')
mye = ad.read_h5ad('../mye_scvi_post_anno.h5ad')
stro = ad.read_h5ad('../Stromal_scVI_post.h5ad')

plot_dotplot(epi,  immune_ligands,   'immune_ligands',       output_dir='../dotplot')
plot_dotplot(lym,  immune_receptors, 'immune_receptors',     output_dir='../dotplot')
plot_dotplot(mye,  immune_ligands,   'immune_ligands_mye',   output_dir='../dotplot')
plot_dotplot(stro, immune_ligands,   'immune_ligands_stro',  output_dir='../dotplot')
