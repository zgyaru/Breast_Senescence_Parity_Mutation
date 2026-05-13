import os
import pandas as pd
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------
# Global settings
# -----------------------------------------------------------------------------
sc.settings.figdir = "../analysis/plots"
os.makedirs(sc.settings.figdir, exist_ok=True)

plt.rcParams.update({
    "font.size": 18,
    "axes.labelsize": 18,
    "xtick.labelsize": 18,
    "ytick.labelsize": 18
})


# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------
def prepare_adata(adata, group_col, order, layer="log1p_norm", exclude=None, rename_map=None):
    """
    Prepare AnnData object for plotting:
    - replace X with a chosen layer
    - rename categories if needed
    - exclude unwanted labels
    - subset to desired order
    - set ordered categorical dtype
    """
    adata = adata.copy()

    if layer is not None:
        if layer not in adata.layers:
            raise KeyError(f"Layer '{layer}' not found in adata.layers")
        adata.X = adata.layers[layer].copy()

    adata.obs[group_col] = adata.obs[group_col].astype(str)

    if rename_map is not None:
        adata.obs[group_col] = adata.obs[group_col].replace(rename_map)

    if exclude is not None:
        adata = adata[~adata.obs[group_col].isin(exclude)].copy()

    adata = adata[adata.obs[group_col].isin(order)].copy()
    adata.obs[group_col] = pd.Categorical(
        adata.obs[group_col],
        categories=order,
        ordered=True
    )

    return adata


def flatten_marker_dict(marker_dict):
    """Flatten dict of marker lists into a unique gene list preserving order."""
    seen = set()
    genes = []
    for gene_list in marker_dict.values():
        for gene in gene_list:
            if gene not in seen:
                seen.add(gene)
                genes.append(gene)
    return genes


def report_missing_genes(adata, marker_dict, label):
    """Print missing genes for quick QC."""
    genes = flatten_marker_dict(marker_dict)
    missing = [g for g in genes if g not in adata.var_names]
    if missing:
        print(f"[{label}] Missing genes ({len(missing)}): {missing}")
    else:
        print(f"[{label}] All marker genes found.")


def make_dotplot(
    adata,
    marker_dict,
    groupby,
    save_scanpy_name=None,
    save_manual_path=None,
    add_totals=False
):
    """
    Generate Scanpy dotplot and save both via scanpy 'save' and manual figure save.
    """
    if save_scanpy_name is not None:
        sc.pl.dotplot(
            adata,
            var_names=marker_dict,
            groupby=groupby,
            dendrogram=False,
            standard_scale="var",
            use_raw=False,
            show=False,
            save=save_scanpy_name
        )

    dp = sc.pl.dotplot(
        adata,
        var_names=marker_dict,
        groupby=groupby,
        dendrogram=False,
        standard_scale="var",
        return_fig=True,
        use_raw=False
    )

    if add_totals:
        dp = dp.add_totals()

    dp.style(dot_edge_color="black", dot_edge_lw=0.5)

    if save_manual_path is not None:
        os.makedirs(os.path.dirname(save_manual_path), exist_ok=True)
        dp.savefig(save_manual_path, dpi=300)


# -----------------------------------------------------------------------------
# Epithelial compartment
# -----------------------------------------------------------------------------
epi_order = [
    "LASP1", "LASP2", "LASP3", "LASP4", "LASP5", "LASP6", "LASP7",
    "LHS1", "LHS2", "LHS3", "LHS4",
    "BMYO1", "BMYO2", "BMYO3", "BMYO4",
    "Tumour", "Aged_Il33+", "DDC"
]

Epithelial_dict = {
    "Epithelial": ["Krt8", "Krt18", "Epcam"],
    "LASP": ["Kit", "Fcgbp", "Lurap1l"],
    "LASP1": ["Fos", "Atf3"],
    "LASP2": ["Csn2", "Csn1s1"],
    "LASP3": ["Csn3", "Cck"],
    "LASP4": ["Aldh1a3", "Ngf"],
    "LASP5": ["Cwh43", "Alox12e"],
    "LASP6": ["Stmn1", "Mki67"],
    "LASP7": ["Ltf", "Ntng1"],
    "LHS": ["Pgr", "Cited1", "Esr1"],
    "LHS1": ["Wfdc18", "Lcn2"],
    "LHS2": ["Prlr", "Nrxn3"],
    "LHS3": ["Krt6a", "Tmem86a"],
    "LHS4": ["Stmn1", "Mki67"],
    "BMYO": ["Trp63", "Krt5", "Krt14"],
    "BMYO1": ["Jag1", "Srm"],
    "BMYO2": ["Egr1", "Gadd45b"],
    "BMYO3": ["Brinp1", "Ptprt"],
    "BMYO4": ["Stmn1", "Mki67"],
    "Tumour": ["Gpx2", "Spp1"],
    "Aged_Il33+": ["Il33"],
    "DDC": ["Erbb4", "Snorc"]
}

epi = ad.read_h5ad("../epi/adata/epi_post_scvi.h5ad")
epi = prepare_adata(
    epi,
    group_col="ct_level3",
    order=epi_order,
    layer="log1p_norm",
    exclude=["Doublet", "Low_Quality"]
)

report_missing_genes(epi, Epithelial_dict, "Epithelial")

make_dotplot(
    epi,
    marker_dict=Epithelial_dict,
    groupby="ct_level3",
    save_scanpy_name="_dotplot_epi.pdf",
    save_manual_path="../analysis/plots/markers_epi.pdf"
)


# -----------------------------------------------------------------------------
# Lymphoid compartment
# -----------------------------------------------------------------------------
lym_order = [
    "T_naive", "CD8_Tcm", "CD8_Tem", "CD8_Trm", "CD8_CTL",
    "CD8_Exhausted", "CD4_Treg", "CD4_Effector", "CD4_Th2",
    "ILC2", "NK", "CD8_NKT", "CD8_Trm/NKT", "DN_NKT",
    "NKT17", "Proliferating_T", "B_cells"
]

lym_dict = {
    "Lymphoid": ["Hcst", "Cd2"],
    "CD3": ["Cd3e", "Cd3d"],
    "CD8": ["Cd8a", "Cd8b1"],
    "Naive/Tcm": ["Sell", "Ccr7"],
    "CD8_Tcm": ["Il7r", "Tcf7", "Cd44"],
    "CD8_Tem": ["Cxcr3", "Gzmk"],
    "CD8_Trm": ["Itgae", "Cxcr6"],
    "CD8_CTL": ["Gzmb", "Tbx21"],
    "CD8_Exhausted": ["Lag3", "Tigit"],
    "CD4": ["Cd4", "Cd40lg"],
    "CD4_Treg": ["Foxp3", "Ctla4"],
    "CD4_Effector": ["Ifng", "Tbx21"],
    "CD4_Th2/ILC2": ["Gata3", "Il5"],
    "NK": ["Ncr1", "Klrk1"],
    "NKT": ["Zbtb16"],
    "NKT17": ["Il17a", "Il17f"],
    "Proliferating_T": ["Mki67"],
    "B_cells": ["Cd74", "Cd79a"]
}

lym = ad.read_h5ad("../lym/adata/lym_post_scvi.h5ad")
lym = prepare_adata(
    lym,
    group_col="ct_level3",
    order=lym_order,
    layer="log1p_norm",
    exclude=["Doublet"]
)

report_missing_genes(lym, lym_dict, "Lymphoid")

make_dotplot(
    lym,
    marker_dict=lym_dict,
    groupby="ct_level3",
    save_scanpy_name="_dotplot_lym.pdf",
    save_manual_path="../analysis/plots/markers_lym.pdf"
)


# -----------------------------------------------------------------------------
# Myeloid compartment
# -----------------------------------------------------------------------------
mye_order = [
    "Mo1", "Mo2", "Mo3_1", "Mo3_2", "Mo3_3", "Tam_1", "Tam_2", "Tam_3",
    "Classical_Monocyte", "Non_Classical_Monocyte", "DC1", "DC2",
    "mDC", "pDC", "Mast Cell", "Neutrophil"
]

mye_dict = {
    "Myeloid": ["Tyrobp", "H2-Aa"],
    "MHCII": ["H2-Ab1", "H2-DMb1"],
    "Macrophage": ["Cd86", "Adgre1"],
    "Mo1/Mo2": ["Itgam"],
    "Mo1": ["Lyve1"],
    "Mo3": ["Mmp12", "Itgax"],
    "Mo3_1": ["Havcr2", "Vcam1"],
    "Mo3_2": ["Apoe", "C1qa"],
    "Mo3_3": ["Fn1", "Retnla"],
    "Tam": ["Spp1"],
    "Tam_1": ["Vegfa", "Cd274"],
    "Tam_2": ["Gatm", "Nme2"],
    "Tam_3": ["Ly6a", "Fcgr4"],
    "Classical_Monocyte": ["Ly6c2", "Ccr2"],
    "Non_Classical_Monocyte": ["Cx3cr1", "Treml4"],
    "Dendritic": ["Zbtb46", "Flt3"],
    "DC1": ["Irf8", "Itgae"],
    "DC2": ["Itgam"],
    "mDC": ["Fscn1", "Ccl22"],
    "pDC": ["Tcf4", "Siglech"],
    "Mast Cell": ["Cpa3", "Tpsb2"],
    "Neutrophil": ["S100a9", "S100a8"]
}

mye = ad.read_h5ad("../mye/adata/mye_scvi_post.h5ad")
mye = prepare_adata(
    mye,
    group_col="ct_level3",
    order=mye_order,
    layer="log1p_norm",
    exclude=["Doublet"]
)

report_missing_genes(mye, mye_dict, "Myeloid")

make_dotplot(
    mye,
    marker_dict=mye_dict,
    groupby="ct_level3",
    save_scanpy_name="_dotplot_mye.pdf",
    save_manual_path="../analysis/plots/markers_mye.pdf",
    add_totals=True
)


# -----------------------------------------------------------------------------
# Stromal compartment
# -----------------------------------------------------------------------------
stro_order = [
    "Fb1", "Fb2", "Fb3", "Fb4", "Fb5", "Fb6", "Fb7", "Fb8",
    "PV1", "PV2", "PV3", "VEC", "VEA1", "VEA2", "VEAT", "LEC1", "LEC2"
]

STRO_dict = {
    "Fibroblast": ["Col3a1", "Dcn"],
    "Fb1": ["Pi16", "Sema3c"],
    "Fb2": ["Fabp4", "Col5a3"],
    "Fb3": ["Gdf10", "Ltbp4"],
    "Fb4": ["Col6a5", "Car4"],
    "Fb5": ["Enpp2", "Rgcc"],
    "Fb6": ["Apod", "Cp"],
    "Fb7": ["Hmcn2", "Myoc"],
    "Fb8": ["Stmn1", "Tubb5"],
    "PV": ["Notch3", "Myh11"],
    "PV1": ["Procr", "Il6"],
    "PV2": ["Rgs5", "Abcc9"],
    "PV3": ["Tagln", "Acta2"],
    "VE": ["Cldn5", "Pecam1"],
    "VEC": ["Plvap", "Emcn"],
    "VEA": ["Sox17", "Efnb2"],
    "VEA1": ["Aqp1", "Rgcc"],
    "VEA2": ["Clu", "Fbln5"],
    "VEAT": ["Pxdn", "Angpt2"],
    "VEV": ["Ackr1", "Selp"],
    "LEC": ["Mmrn1", "Pdpn"],
    "LEC1": ["Plat", "Cfh"],
    "LEC2": ["Lyve1", "Ccl21a"]
}

stro = ad.read_h5ad("../stro/adata/stro_post_scvi.h5ad")
stro = prepare_adata(
    stro,
    group_col="ct_level3",
    order=stro_order,
    layer="log1p_norm",
    exclude=["Doublet"],
    rename_map={
        "VEA_1": "VEA1",
        "VEA_2": "VEA2",
        "LEC_1": "LEC1",
        "LEC_2": "LEC2"
    }
)

report_missing_genes(stro, STRO_dict, "Stromal")

make_dotplot(
    stro,
    marker_dict=STRO_dict,
    groupby="ct_level3",
    save_scanpy_name="_dotplot_stro.pdf",
    save_manual_path="../analysis/plots/markers_stro.pdf"
)
