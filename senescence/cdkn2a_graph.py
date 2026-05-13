import os
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt

from statannotations.Annotator import Annotator
from scipy.stats import mannwhitneyu, linregress, spearmanr
from statsmodels.stats.multitest import multipletests


# =========================
# Paths
# =========================

epi_path = '../new_data/epi_scvi_post.h5ad'

plot_dir_mean = '../sen/plot/new'
plot_dir_prop = '../sen/plot/prop'
plot_dir_main = '../sen/plot'
csv_dir = '../sen/plot/Cdkn2a_proportions'
meta_path = '../analysis/new_senescence/Mouse_ID_Group.csv'

os.makedirs(plot_dir_mean, exist_ok=True)
os.makedirs(plot_dir_prop, exist_ok=True)
os.makedirs(plot_dir_main, exist_ok=True)
os.makedirs(csv_dir, exist_ok=True)


# =========================
# Plot palette / group order
# =========================

palette = {
    'Young_WT_NP': '#000001',
    'Old_WT_NP': '#d28cbc',
    'Old_WT_Parous': '#6acbdd',
    'BRCA1': '#ed1f24',
    'BRCA2': '#008c45',
    'PALB2': '#f7901e',
    'BRCA1_Tumour': '#7f00ff',
    'BRCA2_Tumour': '#ff1493',
    'ABT737_control': '#1f77b4',
    'Treated': '#ffd700'
}

groups = ['Young_WT_NP', 'Old_WT_NP', 'Old_WT_Parous', 'BRCA1', 'BRCA2', 'PALB2']


# =========================
# Helper functions
# =========================

def calculate_prop(adata, gene, expression_threshold=0, group_column='Group', ct_level='ct_level3', sample_col='Mouse_ID'):
    """
    Calculate the proportion of cells expressing a given gene (above threshold)
    grouped by sample and cell type.
    """
    if gene not in adata.var_names:
        raise ValueError(
            f"Gene '{gene}' not found in adata.var_names. Available names include:\n"
            f"{[g for g in adata.var_names if gene.lower() in g.lower()]}"
        )

    gene_expr = adata[:, gene].X
    if not isinstance(gene_expr, np.ndarray):
        gene_expr = gene_expr.toarray()
    gene_expr = gene_expr.flatten()

    gene_col = f'{gene}_expressed'
    obs = adata.obs[[sample_col, ct_level, group_column]].copy()
    obs[gene_col] = (gene_expr > expression_threshold).astype(int)

    df = (
        obs.groupby([sample_col, ct_level, group_column], observed=True)[gene_col]
        .mean()
        .reset_index()
    )
    df.rename(columns={gene_col: f'{gene}_proportion'}, inplace=True)
    df = df.dropna(subset=[f'{gene}_proportion'])
    return df


def cal_mean_expression(adata, gene, group_column='Group', ct_level='ct_level3', sample_col='Mouse_ID'):
    """
    Calculate average expression of `gene` in each sample x cell-type.
    Missing CTs per sample are filled with 0.
    """
    if gene not in adata.var_names:
        close = [g for g in adata.var_names if gene.lower() in g.lower()]
        raise ValueError(
            f"Gene '{gene}' not found in adata.var_names. Close matches: {close[:10]}"
        )

    gx = adata[:, gene].X
    if not isinstance(gx, np.ndarray):
        gx = gx.toarray()
    gx = gx.ravel()

    obs = adata.obs[[sample_col, ct_level, group_column]].copy()
    obs[gene] = gx

    stats = (
        obs.groupby([sample_col, ct_level], observed=True)
        .agg(**{
            f'{gene}_mean_expression': (gene, 'mean'),
            'n_cells': (gene, 'size')
        })
        .reset_index()
    )
    stats['has_cells'] = stats['n_cells'] > 0

    all_samples = obs[sample_col].unique()
    all_cts = sorted(obs[ct_level].dropna().unique())
    full = pd.MultiIndex.from_product(
        [all_samples, all_cts],
        names=[sample_col, ct_level]
    ).to_frame(index=False)

    out = full.merge(stats, on=[sample_col, ct_level], how='left')

    sample_to_group = (
        obs[[sample_col, group_column]]
        .drop_duplicates(subset=[sample_col])
        .set_index(sample_col)[group_column]
    )
    out[group_column] = out[sample_col].map(sample_to_group)

    out[f'{gene}_mean_expression'] = out[f'{gene}_mean_expression'].fillna(0.0)
    out['n_cells'] = out['n_cells'].fillna(0).astype(int)
    out['has_cells'] = out['has_cells'].fillna(False)

    cols = [
        sample_col, ct_level, group_column,
        f'{gene}_mean_expression', 'n_cells', 'has_cells'
    ]
    return out[cols].sort_values([group_column, sample_col, ct_level]).reset_index(drop=True)


def perform_test(df, group_column='Group', gene='Cdkn2a', alpha=0.05, adjust='holm',
                 ref_group='Young_WT_NP', alternative='two-sided', dropna=True):
    
    if gene not in df.columns or group_column not in df.columns:
        raise ValueError(f"Columns '{group_column}' and/or '{gene}' not found.")

    sub = df[[group_column, gene]].copy()
    sub[gene] = pd.to_numeric(sub[gene], errors='coerce')
    if dropna:
        sub = sub.dropna(subset=[group_column, gene])

    if ref_group not in set(sub[group_column].astype(str)):
        raise ValueError(f"Reference group '{ref_group}' not found in '{group_column}'.")

    groups_dict = {g: v[gene].values for g, v in sub.groupby(group_column)}
    ref_vals = groups_dict.get(ref_group)

    if ref_vals is None or len(ref_vals) == 0:
        raise ValueError(f"No data for reference group '{ref_group}'.")

    comparisons = []
    for g, vals in groups_dict.items():
        if g == ref_group:
            continue
        if len(vals) == 0 or len(ref_vals) == 0:
            continue

        U, p = mannwhitneyu(ref_vals, vals, alternative=alternative)
        n1, n2 = len(ref_vals), len(vals)
        r_rb = 1 - (2 * U) / (n1 * n2)

        comparisons.append({
            'group1': ref_group,
            'group2': g,
            'n1': n1,
            'n2': n2,
            'U': U,
            'p_raw': p,
            'effect_rbs': r_rb
        })

    if not comparisons:
        return pd.DataFrame(columns=['group1', 'group2', 'n1', 'n2', 'U', 'p_raw', 'p_adj', 'reject', 'effect_rbs'])

    out = pd.DataFrame(comparisons)

    if adjust is None:
        out['p_adj'] = out['p_raw']
        out['reject'] = out['p_adj'] < alpha
    else:
        if adjust not in {'holm', 'bonferroni', 'fdr_bh'}:
            raise ValueError("adjust must be one of {'holm','bonferroni','fdr_bh', None}")
        rej, p_adj, _, _ = multipletests(out['p_raw'].values, method=adjust, alpha=alpha)
        out['p_adj'] = p_adj
        out['reject'] = rej

    out = out.sort_values('p_adj').reset_index(drop=True)
    return out


def plot_gene_proportions(df, gene, order=None, group_column='Group', ct_level='ct_level3',
                          palette=None, save=None, data=None):
    if order is None:
        order = sorted(df[group_column].dropna().unique())

    sns.set(style="whitegrid", context="talk", font_scale=1.1)

    if data == 'mean':
        value_col = f"{gene}_mean_expression"
    else:
        value_col = f"{gene}_proportion"

    unique_cts = sorted(df[ct_level].dropna().unique())
    num_groups = len(order)

    gene_dir = os.path.join(save, f'{gene}_proportions')
    os.makedirs(gene_dir, exist_ok=True)

    for ct in unique_cts:
        subset = df[df[ct_level] == ct].copy()
        if subset.empty:
            continue

        fig_width = max(1 * num_groups, 5)
        plt.figure(figsize=(fig_width, 6))

        ax = sns.barplot(
            data=subset,
            x=group_column,
            y=value_col,
            errorbar='sd',
            capsize=0.1,
            palette=palette,
            errcolor='black',
            errwidth=1.2,
            order=order
        )

        sns.swarmplot(
            data=subset,
            x=group_column,
            y=value_col,
            color='black',
            size=3,
            alpha=0.8,
            dodge=True,
            order=order
        )

        pairs = [
            ('Young_WT_NP', 'Old_WT_NP'),
            ('Young_WT_NP', 'Old_WT_Parous'),
            ('Young_WT_NP', 'BRCA1'),
            ('Young_WT_NP', 'BRCA2'),
            ('Young_WT_NP', 'PALB2')
        ]
        present_groups = set(subset[group_column].dropna().unique())
        pairs = [p for p in pairs if p[0] in present_groups and p[1] in present_groups]

        if pairs:
            annotator = Annotator(ax, pairs, data=subset, x=group_column, y=value_col, order=order)
            annotator.configure(
                test='Mann-Whitney',
                text_format='star',
                loc='inside',
                verbose=0,
                comparisons_correction='BH'
            )
            annotator.apply_and_annotate()

        ax.set_ylabel(f'{gene}+ {data}', fontsize=12)
        ax.set_xlabel(group_column, fontsize=12)
        ax.tick_params(axis='x', labelrotation=90)

        ymax = subset[value_col].max()
        upper = max(0.6, ymax * 1.25 if pd.notna(ymax) else 0.6)
        ax.set_ylim(0, upper)

        sns.despine()
        plt.tight_layout()

        ct_str = str(ct).replace('/', '_').replace(' ', '_')
        save_path_pdf = os.path.join(gene_dir, f'{gene}_{ct_str}_{group_column}_{data}.pdf')
        save_path_svg = os.path.join(gene_dir, f'{gene}_{ct_str}_{group_column}_{data}.svg')

        print(f"Saving figure to {save_path_pdf}")
        plt.savefig(save_path_pdf, dpi=300, transparent=True)
        plt.savefig(save_path_svg, dpi=300, transparent=True)
        plt.close()


def plot_scatter_with_regression(
    df,
    x_col,
    y_col,
    title,
    save_path,
    pos=(0.05, 0.95),
    decimals=3,
    ci=None,
    spearman=False,
    y_is_percentage=False,
    auto_scale_proportion=True,
    clip_to_bounds=True,
    y_min=0.0,
    y_max=100.0,
    note_on_plot=True,
    dpi=600
):
    d = df[[x_col, y_col]].replace([np.inf, -np.inf], np.nan).dropna()
    if len(d) < 2:
        raise ValueError("Not enough valid points after cleaning (need ≥2).")

    d = d.copy()
    scaled = False
    clipped = False
    y_label = y_col

    if y_is_percentage:
        y_vals = d[y_col].to_numpy(dtype=float)
        if auto_scale_proportion and np.nanmin(y_vals) >= 0 and np.nanmax(y_vals) <= 1.0:
            d[y_col] = y_vals * 100.0
            scaled = True

        if clip_to_bounds:
            before = d[y_col].to_numpy(dtype=float)
            after = np.clip(before, y_min, y_max)
            if np.any(after != before):
                clipped = True
            d[y_col] = after

        if "%" not in y_label:
            y_label = f"{y_label} (%)"

        if not clip_to_bounds:
            ymin_observed = float(np.nanmin(d[y_col]))
            ymax_observed = float(np.nanmax(d[y_col]))
            if ymin_observed < y_min - 1e-12 or ymax_observed > y_max + 1e-12:
                raise ValueError(
                    f"{y_col} contains values outside [{y_min}, {y_max}] but clip_to_bounds=False. "
                    f"Observed range: [{ymin_observed:.3g}, {ymax_observed:.3g}]."
                )

    x = d[x_col].to_numpy(dtype=float)
    y = d[y_col].to_numpy(dtype=float)

    if np.allclose(np.nanstd(x), 0):
        raise ValueError("x is constant; regression is undefined.")
    if np.allclose(np.nanstd(y), 0):
        raise ValueError("y is constant; regression is undefined.")

    slope, intercept, r_value, p_value, std_err = linregress(x, y)
    R2 = r_value ** 2

    spearman_str = ""
    if spearman:
        rho, p_s = spearmanr(x, y)
        spearman_str = f"\n$\\rho_s$ = {rho:.{decimals}f} (P = {p_s:.2e})"

    annotation = (
        fr"$\beta$ = {slope:.{decimals}g}"
        f"\n$R^2$ = {R2:.{decimals}f}"
        f"\nP = {p_value:.2e}"
        f"{spearman_str}"
    )

    sns.set_style("white")
    sns.set_context("paper", font_scale=1.6)
    plt.figure(figsize=(6, 5))

    ax = sns.regplot(
        x=x_col, y=y_col, data=d,
        scatter_kws={'s': 40, 'alpha': 0.7, 'color': 'black', 'edgecolor': 'none'},
        line_kws={'linewidth': 2, 'color': 'black'},
        ci=ci
    )

    for spine in ('top', 'right'):
        ax.spines[spine].set_visible(False)
    for spine in ('left', 'bottom'):
        ax.spines[spine].set_linewidth(1)

    ax.set_xlabel(x_col, fontsize=14)
    ax.set_ylabel(y_label if y_is_percentage else y_col, fontsize=14)
    ax.set_title(title, fontsize=16, pad=10)

    ax.text(
        pos[0], pos[1], annotation, transform=ax.transAxes,
        fontsize=12, va='top', ha='left',
        bbox=dict(boxstyle='round,pad=0.25', fc='white', ec='none', alpha=0.9)
    )

    if note_on_plot and (scaled or clipped):
        note_lines = []
        if scaled:
            note_lines.append("Scaled: proportion → % (×100)")
        if clipped:
            note_lines.append(f"Clipped to [{y_min:g}, {y_max:g}]")
        ax.text(
            0.99, 0.01, " | ".join(note_lines),
            transform=ax.transAxes, ha='right', va='bottom',
            fontsize=9, color='gray'
        )

    plt.tight_layout()
    plt.savefig(save_path, dpi=dpi, bbox_inches='tight', transparent=True)
    plt.close()

    stats = {
        "slope": slope,
        "intercept": intercept,
        "r": r_value,
        "R2": R2,
        "p": p_value,
        "std_err": std_err,
        "scaled_proportion_to_percent": scaled if y_is_percentage else False,
        "clipped_to_bounds": clipped if y_is_percentage else False
    }
    if spearman:
        stats.update({"spearman_rho": rho, "spearman_p": p_s})

    return stats


# =========================
# Load / filter epithelial data
# =========================

epi = ad.read_h5ad(epi_path)
epi.obs['Group'] = epi.obs['Group'].astype(str)

epi = epi[epi.obs['Genotype'].isin(['WT', 'BRCA1', 'PALB2', 'BRCA2'])].copy()
epi = epi[~epi.obs['Mouse_ID'].isin(['SIGAH2'])].copy()
epi = epi[epi.obs['ct_level2'].isin(['LASP'])].copy()
epi = epi[epi.obs['ct_level3'] != 'Doublet'].copy()
epi = epi[epi.obs['ct_level3'] != 'DDC'].copy()
epi = epi[epi.obs['ct_level3'] != 'Low_Quality'].copy()
epi = epi[epi.obs['Parity'].isin(['NP', 'Parous'])].copy()

if pd.api.types.is_categorical_dtype(epi.obs['Mouse_ID']):
    epi.obs['Mouse_ID'] = epi.obs['Mouse_ID'].cat.remove_unused_categories()


# =========================
# Mean expression plots
# =========================

df_mean_level2 = cal_mean_expression(epi, 'Cdkn2a', group_column='Group', ct_level='ct_level2')
df_mean_level2 = df_mean_level2[df_mean_level2['Group'].isin(groups)].copy()

plot_gene_proportions(
    df_mean_level2,
    'Cdkn2a',
    order=groups,
    group_column='Group',
    ct_level='ct_level2',
    palette=palette,
    save=plot_dir_mean,
    data='mean'
)
df_mean_level2.to_csv(os.path.join(csv_dir, 'Cdkn2a_mean_level2.csv'), index=False)

df_mean_level3 = cal_mean_expression(epi, 'Cdkn2a', group_column='Group', ct_level='ct_level3')
df_mean_level3 = df_mean_level3[df_mean_level3['Group'].isin(groups)].copy()

plot_gene_proportions(
    df_mean_level3,
    'Cdkn2a',
    order=groups,
    group_column='Group',
    ct_level='ct_level3',
    palette=palette,
    save=plot_dir_mean,
    data='mean'
)
df_mean_level3.to_csv(os.path.join(csv_dir, 'Cdkn2a_mean_level3.csv'), index=False)


# =========================
# Proportion plots
# =========================

df_prop_level2 = calculate_prop(epi, 'Cdkn2a', expression_threshold=0, group_column='Group', ct_level='ct_level2')
df_prop_level2 = df_prop_level2[df_prop_level2['Group'].isin(groups)].copy()

plot_gene_proportions(
    df_prop_level2,
    'Cdkn2a',
    order=groups,
    group_column='Group',
    ct_level='ct_level2',
    palette=palette,
    save=plot_dir_prop,
    data='proportion'
)

df_prop_level2.to_csv(os.path.join(csv_dir, 'Cdkn2a_proportions_level2.csv'), index=False)


df_prop_level3 = calculate_prop(epi, 'Cdkn2a', expression_threshold=0, group_column='Group', ct_level='ct_level3')
df_prop_level3 = df_prop_level3[df_prop_level3['Group'].isin(groups)].copy()
plot_gene_proportions(
    df_prop_level3,
    'Cdkn2a',
    order=groups,
    group_column='Group',
    ct_level='ct_level3',
    palette=palette,
    save=plot_dir_prop,
    data='proportion'
)

df_prop_level3.to_csv(os.path.join(csv_dir, 'Cdkn2a_proportions_level3.csv'), index=False)

# =========================
# Statistical tests
# =========================

result_mean_level2 = perform_test(
    df_mean_level2,
    group_column='Group',
    gene='Cdkn2a_mean_expression',
    alpha=0.05,
    adjust='fdr_bh',
    ref_group='Young_WT_NP',
    alternative='two-sided',
    dropna=True
)
result_mean_level2.to_csv(os.path.join(csv_dir, 'Cdkn2a_mean_level2_test.csv'), index=False)

LASP2 = df_mean_level3[df_mean_level3['ct_level3'] == 'LASP2'].copy()
lasp2_result = perform_test(
    LASP2,
    group_column='Group',
    gene='Cdkn2a_mean_expression',
    alpha=0.05,
    adjust='fdr_bh',
    ref_group='Young_WT_NP',
    alternative='two-sided',
    dropna=True
)
lasp2_result.to_csv(os.path.join(csv_dir, 'Cdkn2a_mean_LASP2_test.csv'), index=False)

LASP5 = df_mean_level3[df_mean_level3['ct_level3'] == 'LASP5'].copy()
lasp5_result = perform_test(
    LASP5,
    group_column='Group',
    gene='Cdkn2a_mean_expression',
    alpha=0.05,
    adjust='fdr_bh',
    ref_group='Young_WT_NP',
    alternative='two-sided',
    dropna=True
)
lasp5_result.to_csv(os.path.join(csv_dir, 'Cdkn2a_mean_LASP5_test.csv'), index=False)


# =========================
# Age metadata merge for proportion regression
# =========================

proportion_df = calculate_prop(epi, 'Cdkn2a', expression_threshold=0, group_column='Group', ct_level='ct_level2')
proportion_df = proportion_df.reset_index(drop=True)

meta = pd.read_csv(meta_path)
meta = meta[['Mouse_ID', 'Age']]

proportion_df = proportion_df.merge(meta, on='Mouse_ID', how='left')

brca1 = proportion_df[proportion_df['Group'] == 'BRCA1'].copy()
brca2 = proportion_df[proportion_df['Group'] == 'BRCA2'].copy()
palb2 = proportion_df[proportion_df['Group'] == 'PALB2'].copy()
wt = proportion_df[proportion_df['Group'].isin(['Young_WT_NP', 'Old_WT_NP'])].copy()


# =========================
# Proportion vs age regression
# =========================

plot_scatter_with_regression(
    brca1, 'Age', 'Cdkn2a_proportion', 'Cdkn2a Proportion in BRCA1',
    os.path.join(plot_dir_main, 'Cdkn2a_Brca1_prop.svg')
)
plot_scatter_with_regression(
    brca2, 'Age', 'Cdkn2a_proportion', 'Cdkn2a Proportion in BRCA2',
    os.path.join(plot_dir_main, 'Cdkn2a_Brca2_prop.svg')
)
plot_scatter_with_regression(
    palb2, 'Age', 'Cdkn2a_proportion', 'Cdkn2a Proportion in PALB2',
    os.path.join(plot_dir_main, 'Cdkn2a_Palb2_prop.svg')
)
plot_scatter_with_regression(
    wt, 'Age', 'Cdkn2a_proportion', 'Cdkn2a Proportion in WT',
    os.path.join(plot_dir_main, 'Cdkn2a_WT_prop.svg')
)

plot_scatter_with_regression(
    brca1, 'Age', 'Cdkn2a_proportion', 'Cdkn2a Proportion in BRCA1',
    os.path.join(plot_dir_main, 'Cdkn2a_Brca1_prop.pdf')
)
plot_scatter_with_regression(
    brca2, 'Age', 'Cdkn2a_proportion', 'Cdkn2a Proportion in BRCA2',
    os.path.join(plot_dir_main, 'Cdkn2a_Brca2_prop.pdf')
)
plot_scatter_with_regression(
    palb2, 'Age', 'Cdkn2a_proportion', 'Cdkn2a Proportion in PALB2',
    os.path.join(plot_dir_main, 'Cdkn2a_Palb2_prop.pdf')
)
plot_scatter_with_regression(
    wt, 'Age', 'Cdkn2a_proportion', 'Cdkn2a Proportion in WT',
    os.path.join(plot_dir_main, 'Cdkn2a_WT_prop.pdf')
)


# =========================
# Mean expression vs age regression
# =========================

level2 = pd.read_csv(os.path.join(csv_dir, 'Cdkn2a_mean_level2.csv'))
meta = pd.read_csv(meta_path)
meta = meta[['Mouse_ID', 'Age']]

level2 = level2.merge(meta, on='Mouse_ID', how='left')

LASP_Brca1 = level2[level2['Group'] == 'BRCA1'].copy()
LASP_Brca2 = level2[level2['Group'] == 'BRCA2'].copy()
LASP_Palb2 = level2[level2['Group'] == 'PALB2'].copy()
LASP_WT = level2[level2['Group'].isin(['Young_WT_NP', 'Old_WT_NP'])].copy()

plot_scatter_with_regression(
    LASP_Brca1, 'Age', 'Cdkn2a_mean_expression',
    'Cdkn2a Proportion in LASP BRCA1',
    os.path.join(plot_dir_main, 'Cdkn2a_LASP_Brca1_new.pdf')
)
plot_scatter_with_regression(
    LASP_Brca2, 'Age', 'Cdkn2a_mean_expression',
    'Cdkn2a Proportion in LASP BRCA2',
    os.path.join(plot_dir_main, 'Cdkn2a_LASP_Brca2_new.pdf')
)
plot_scatter_with_regression(
    LASP_Palb2, 'Age', 'Cdkn2a_mean_expression',
    'Cdkn2a Proportion in LASP PALB2',
    os.path.join(plot_dir_main, 'Cdkn2a_LASP_Palb2_new.pdf')
)
plot_scatter_with_regression(
    LASP_WT, 'Age', 'Cdkn2a_mean_expression',
    'Cdkn2a Proportion in LASP WT',
    os.path.join(plot_dir_main, 'Cdkn2a_LASP_WT_new.pdf')
)

plot_scatter_with_regression(
    LASP_Brca1, 'Age', 'Cdkn2a_mean_expression',
    'Cdkn2a Proportion in LASP BRCA1',
    os.path.join(plot_dir_main, 'Cdkn2a_LASP_Brca1_mean.svg')
)
plot_scatter_with_regression(
    LASP_Brca2, 'Age', 'Cdkn2a_mean_expression',
    'Cdkn2a Proportion in LASP BRCA2',
    os.path.join(plot_dir_main, 'Cdkn2a_LASP_Brca2_mean.svg')
)
plot_scatter_with_regression(
    LASP_Palb2, 'Age', 'Cdkn2a_mean_expression',
    'Cdkn2a Proportion in LASP PALB2',
    os.path.join(plot_dir_main, 'Cdkn2a_LASP_Palb2_mean.svg')
)
plot_scatter_with_regression(
    LASP_WT, 'Age', 'Cdkn2a_mean_expression',
    'Cdkn2a Proportion in LASP WT',
    os.path.join(plot_dir_main, 'Cdkn2a_LASP_WT_mean.svg')
)

