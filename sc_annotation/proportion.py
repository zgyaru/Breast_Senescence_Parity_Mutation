import os
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns

figdir = '../marker_cluster_dotplot/figures/proportion'
os.makedirs(figdir, exist_ok=True)


# Nature-like style settings
mpl.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 300,
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "axes.linewidth": 0.6,
    "pdf.fonttype": 42,
    "ps.fonttype": 42
})


def extract_proportions(adata, groupby='Group', ct_level='ct_level3', ct_order=None, group_order=None):
    obs = adata.obs.copy()

    if ct_order is None:
        ct_order = obs[ct_level].value_counts().index.tolist()
    if group_order is None:
        group_order = obs[groupby].dropna().unique().tolist()

    obs = obs.dropna(subset=[groupby, ct_level]).copy()
    obs[ct_level] = pd.Categorical(obs[ct_level], categories=ct_order, ordered=True)
    obs[groupby] = pd.Categorical(obs[groupby], categories=group_order, ordered=True)

    counts = (
        obs.groupby([groupby, ct_level], observed=False)
        .size()
        .unstack(fill_value=0)
        .reindex(index=group_order, columns=ct_order, fill_value=0)
    )
    props = counts.div(counts.sum(axis=1), axis=0).fillna(0)

    return props, counts


def plot_stacked_bars(proportions, width_mm=85, height_mm=55, save_path=None, title=None):
    groups = proportions.index.tolist()
    cts = proportions.columns.tolist()

    base_palette = sns.color_palette("tab20", 20)
    palette_map = {ct: base_palette[i % len(base_palette)] for i, ct in enumerate(cts)}

    fig_w, fig_h = width_mm / 25.4, height_mm / 25.4
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    x = np.arange(len(groups))
    bottoms = np.zeros(len(groups))

    for ct in cts:
        vals = proportions[ct].values
        ax.bar(
            x, vals, bottom=bottoms, width=0.8,
            color=palette_map[ct], label=ct, edgecolor="none"
        )
        bottoms += vals

    ax.set_ylim(0, 1)
    ax.set_ylabel("Proportion")
    ax.set_xticks(x)
    ax.set_xticklabels(groups, rotation=45, ha="right")
    ax.set_yticks(np.arange(0, 1.01, 0.2))

    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)

    ax.grid(axis="y", linewidth=0.3, alpha=0.5)
    ax.tick_params(axis="both", length=3, width=0.6)

    legend = ax.legend(
        title='Cell Type',
        bbox_to_anchor=(1.05, 1),
        loc="upper left",
        frameon=False,
        handlelength=1.2,
        columnspacing=0.8
    )

    if title:
        ax.set_title(title, pad=6)

    fig.tight_layout(pad=0.8)

    if save_path:
        base, _ = os.path.splitext(save_path)
        fig.savefig(base + ".svg", bbox_inches="tight", bbox_extra_artists=[legend])
        fig.savefig(base + ".png", dpi=300, bbox_inches="tight", bbox_extra_artists=[legend])

    return fig, ax


# =========================
# Preprocess and save
# =========================

epi = ad.read_h5ad('../epi/adata/epi_scvi_post.h5ad')
epi.obs['ct_level3'] = epi.obs['ct_level3'].astype(str)
epi.obs.loc[epi.obs['ct_level3'] == 'Aged_Il33+', 'ct_level3'] = 'Il33+'
epi.write_h5ad('../epi/adata/epi_scvi_post.h5ad')

mye = ad.read_h5ad('../mye/adata/mye_scvi_post.h5ad')
mye.obs['ct_level3'] = mye.obs['ct_level3'].astype(str)
mye.obs.loc[mye.obs['ct_level3'] == 'Tam_1', 'ct_level3'] = 'Tam1'
mye.obs.loc[mye.obs['ct_level3'] == 'Tam_2', 'ct_level3'] = 'Tam2'
mye.obs.loc[mye.obs['ct_level3'] == 'Tam_3', 'ct_level3'] = 'Tam3'
mye.write_h5ad('../mye/adata/mye_scvi_post_anno.h5ad')

stro = ad.read_h5ad('../stromal/adata/stromal_scvi_post.h5ad')
stro.obs['ct_level3'] = stro.obs['ct_level3'].astype(str)
stro.obs.loc[stro.obs['ct_level3'] == 'Fb1', 'ct_level3'] = 'FB1'
stro.obs.loc[stro.obs['ct_level3'] == 'Fb2', 'ct_level3'] = 'FB2'
stro.obs.loc[stro.obs['ct_level3'] == 'Fb3', 'ct_level3'] = 'FB3'
stro.obs.loc[stro.obs['ct_level3'] == 'Fb4', 'ct_level3'] = 'FB4'
stro.obs.loc[stro.obs['ct_level3'] == 'Fb5', 'ct_level3'] = 'FB5'
stro.obs.loc[stro.obs['ct_level3'] == 'Fb6', 'ct_level3'] = 'FB6'
stro.obs.loc[stro.obs['ct_level3'] == 'Fb7', 'ct_level3'] = 'FB7'
stro.obs.loc[stro.obs['ct_level3'] == 'Fb8', 'ct_level3'] = 'FB8'
stro.obs.loc[stro.obs['ct_level3'] == 'VEA_1', 'ct_level3'] = 'VEA1'
stro.obs.loc[stro.obs['ct_level3'] == 'VEA_2', 'ct_level3'] = 'VEA2'
stro.write_h5ad('../stromal/adata/stromal_scvi_post.h5ad')


# Shared order
group_order = ['Young_WT_NP', 'Old_WT_NP', 'Old_WT_Parous', 'BRCA1', 'BRCA2', 'PALB2']


# =========================
# Epithelial
# =========================

epi = ad.read_h5ad('../epi/adata/epi_scvi_post.h5ad')
epi = epi[epi.obs['Parity'].isin(['NP', 'Parous'])].copy()
epi = epi[~epi.obs['ct_level3'].isin(['Doublet', 'Low_Quality'])].copy()
epi = epi[epi.obs['Group'].isin(group_order)].copy()

epi.obs['plot_level3'] = epi.obs['ct_level3'].astype(str)
epi.obs.loc[epi.obs['ct_level3'].isin(['LHS1', 'LHS2', 'LHS3', 'LHS4']), 'plot_level3'] = 'LHS'
epi.obs.loc[epi.obs['ct_level3'].isin(['BMYO1', 'BMYO2', 'BMYO3', 'BMYO4']), 'plot_level3'] = 'BMYO'

ct_order = ['LASP1', 'LASP2', 'LASP3', 'LASP4', 'LASP5', 'LASP6', 'LASP7', 'LHS', 'BMYO', 'Tumour', 'Il33+', 'DDC']
epi_prop, epi_counts = extract_proportions(
    epi, groupby='Group', ct_level='plot_level3',
    ct_order=ct_order, group_order=group_order
)

plot_stacked_bars(
    epi_prop,
    width_mm=85, height_mm=80,
    save_path=os.path.join(figdir, 'epi_proportion_bar_plot.svg')
)
epi_prop.to_csv(os.path.join(figdir, 'epi_proportion_bar_plot.csv'))


# =========================
# Stromal
# =========================

stro = ad.read_h5ad('../stromal/adata/stromal_scvi_post.h5ad')
stro = stro[stro.obs['Parity'].isin(['NP', 'Parous'])].copy()
stro = stro[~stro.obs['ct_level3'].isin(['Doublet', 'Low_Quality'])].copy()
stro = stro[stro.obs['Group'].isin(group_order)].copy()

stro.obs['plot_level3'] = stro.obs['ct_level3'].astype(str)
stro.obs.loc[~stro.obs['ct_level3'].isin(['FB1', 'FB2', 'FB3', 'FB4', 'FB5', 'FB6', 'FB7', 'FB8']), 'plot_level3'] = 'Vascular'

ct_order = ['FB1', 'FB2', 'FB3', 'FB4', 'FB5', 'FB6', 'FB7', 'FB8', 'Vascular']
stro_prop, stro_counts = extract_proportions(
    stro, groupby='Group', ct_level='plot_level3',
    ct_order=ct_order, group_order=group_order
)

plot_stacked_bars(
    stro_prop,
    width_mm=85, height_mm=80,
    save_path=os.path.join(figdir, 'stro_proportion_bar_plot.svg')
)
stro_prop.to_csv(os.path.join(figdir, 'stro_proportion_bar_plot.csv'))


# =========================
# Myeloid
# =========================

mye = ad.read_h5ad('../mye/adata/mye_scvi_post.h5ad')
mye = mye[mye.obs['Parity'].isin(['NP', 'Parous'])].copy()
mye = mye[~mye.obs['ct_level3'].isin(['Doublet', 'Low_Quality'])].copy()
mye = mye[mye.obs['Group'].isin(group_order)].copy()

ct_order = [
    'Mo1', 'Mo2', 'Mo3_1', 'Mo3_2', 'Mo3_3', 'Tam1', 'Tam2', 'Tam3',
    'Classical_Monocyte', 'Non_Classical_Monocyte', 'DC1', 'DC2',
    'mDC', 'pDC', 'Mast Cell', 'Neutrophil'
]
mye_prop, mye_counts = extract_proportions(
    mye, groupby='Group', ct_level='ct_level3',
    ct_order=ct_order, group_order=group_order
)

plot_stacked_bars(
    mye_prop,
    width_mm=120, height_mm=80,
    save_path=os.path.join(figdir, 'mye_proportion_bar_plot.svg')
)
mye_prop.to_csv(os.path.join(figdir, 'mye_proportion_bar_plot.csv'))


# =========================
# Lymphoid
# =========================

lym = ad.read_h5ad('../lym/adata/lym_post_scvi.h5ad')
lym = lym[lym.obs['Parity'].isin(['NP', 'Parous'])].copy()
lym = lym[~lym.obs['ct_level3'].isin(['Doublet', 'Low_Quality'])].copy()
lym = lym[lym.obs['Group'].isin(group_order)].copy()

ct_order = [
    "T_naive", "CD8_Tcm", "CD8_Tem", "CD8_Trm", "CD8_CTL",
    "CD8_Exhausted", "CD4_Treg", "CD4_Effector", "CD4_Th2", "ILC2",
    "NK", "CD8_NKT", 'CD8_Trm/NKT', "DN_NKT", "NKT17",
    "Proliferating_T", "B_cells"
]
lym_prop, lym_counts = extract_proportions(
    lym, groupby='Group', ct_level='ct_level3',
    ct_order=ct_order, group_order=group_order
)

plot_stacked_bars(
    lym_prop,
    width_mm=120, height_mm=80,
    save_path=os.path.join(figdir, 'lym_proportion_bar_plot_final.svg')
)
lym_prop.to_csv(os.path.join(figdir, 'lym_proportion_bar_plot.csv'))
