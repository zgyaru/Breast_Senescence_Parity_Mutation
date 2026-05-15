import liana as li
import scanpy as sc
import pandas as pd
#import squidpy as sq
import os
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42      # using TrueType
mpl.rcParams['ps.fonttype'] = 42 
mpl.rcParams['font.sans-serif'] = 'Arial'
plt.rcParams["font.monospace"] = 'Arial'
folder = '/mnt/home5/pharmacology/yz2071/scratch/projects/Breast-mouse'
os.chdir(folder)

data = sc.read_h5ad('./data/visium/annotated/WKBR_adata_tacco_anno.h5ad')
br = pd.read_csv('./results/BR/data4muspan_BR_region5.csv',index_col=0)
data.obs['region'] = br.loc[data.obs_names]['region']

li.ut.spatial_neighbors(data, bandwidth=15, cutoff=0.05, kernel='gaussian', set_diag=False)

lrdata = li.mt.bivariate(data,
                    resource_name='mouseconsensus',
                    local_name='morans', # Name of the function
                    global_name="morans", # Name global function
                    n_perms=1000, # Number of permutations to calculate a p-value
                    mask_negatives=False, # Whether to mask LowLow/NegativeNegative interactions
                    add_categories=True, # Whether to add local categories to the results
                    nz_prop=0, # Minimum expr. proportion for ligands/receptors and their subunits
                    use_raw=False,
                    verbose=True
                    )

lrdata.write_h5ad('./results/BR/liana_noMask_np.h5ad')


from liana.method import singlecellsignalr, connectome, cellphonedb, natmi, logfc, cellchat, geometric_mean
regions = data.obs['region'].unique().tolist()
regions = [x for x in regions if str(x) != 'nan']
inter_data = li.rs.select_resource('mouseconsensus')
inter_data = pd.concat([pd.DataFrame([['Vtcn1','Pdcd1']], columns=inter_data.columns), inter_data])

result_list = {}
for i in regions:
    print(i)
    br_sub = data[data.obs['region'] == i].copy()
    cellphonedb(br_sub,
            groupby='Tacco_final',
            #resource_name='mouseconsensus',
            resource = inter_data,
            expr_prop=0.01,
            verbose=True, 
            key_added='cpdb_res',
            use_raw=False)
    br_sub.uns['cpdb_res']['cluster'] = i
    result_list[i] = br_sub.uns['cpdb_res']

result = pd.concat(result_list.values())
result.to_csv('./results/BR/BR_cluster_cellphone_db_addVtcn1.csv')



import numpy as np
#df = result.copy()
df = pd.read_csv('./results/BR/BR_cluster_cellphone_db_addVtcn1.csv')
df['cell_pair'] = df['source'] + ':' + df['target']
df['interaction'] = df['ligand'] + '^' + df['receptor']
df['log10_pval'] = -np.log10(df['cellphone_pvals'].replace(0, 1e-300))
df['log2_mean'] = np.log2(df['lr_means'].replace(0, np.nan))

df['source'].unique().tolist()
df['target'].unique().tolist()
df['cell_pair'].unique().tolist()

pairs = ['Mac/Mono:T_cells', 'Tam:T_cells', 'LASP2:T_cells','LASP2:Tam', 'LASP2:Mac/Mono', 'LASP5:T_cells', 'LASP5:Tam', 'LASP5:Mac/Mono', 
         'Tumour-like:T_cells', 'Tumour-like:Tam', 'Tumour-like:Mac/Mono', 'LASP_Other:T_cells', 'LASP_Other:Tam', 'LASP_Other:Mac/Mono']


pairs = ['Mac/Mono:T_cells', 
         'Tam:T_cells', 
         'LASP2:T_cells',
         'LASP2:Tam', 
         'LASP2:Mac/Mono', 
         'LASP2:DC', 
         'LASP2:B_cells',
         'LASP2:Fibroblasts',
         'LASP5:T_cells', 
         'LASP5:Tam', 
         'LASP5:Mac/Mono', 
         'LASP5:B_cells', 
         'LASP5:DC', 
         'LASP5:Fibroblasts',
         'Tumour-like:T_cells', 
         'Tumour-like:Tam', 
         'Tumour-like:Mac/Mono', 
         'Tumour-like:B_cells', 
         'Tumour-like:DC',
         'Tumour-like:Fibroblasts',
         'LASP_Other:T_cells', 
         'LASP_Other:Tam', 
         'LASP_Other:Mac/Mono',
         'LASP_Other:B_cells', 
         'LASP_Other:DC', 
         'LASP_Other:Fibroblasts',
         'Fibroblast:Tumour-like',
         'Fibroblast:LASP2',
         'Fibroblast:LASP5',
         'Fibroblast:LASP_Other'
         ]


region_top_pairs = [
    "Mfge8^Itgb3",   ## R1
    "App^Ptger2",
    "C3^C5ar2",
    "C4b^C5ar2",
    "Cd40lg^Cd9",
    "Wnt5b^Lrp6",
    "H2-K1^Erbb2",   ## R2
    "Wnt5b^Fzd4",
    "Col4a1^Itga2",
    "Col1a2^Itga2b",
    "Il34^Ptprz1",
    "Ccl5^Ccr1",
    "App^Ngfr",    ## R3
    "Spp1^Itga4",
    "Thbs2^Cd47",
    "Mfge8^Itgb3",
    "Fn1^Itga4",
    "Ccl5^Ccrl2",
    "Ltf^Tnfrsf11b",    ## R4
    "App^Lrp6",
    "Ltf^Tnfrsf11b",
    "Col1a2^Itgb8",
    "Col1a2^Itga1",
    "Tnc^Sdc4",
    "Lgals3^Eng",    ## R5
    "Col3a1^Itga1",
    "Col1a1^Itga1",
    "Wnt5b^Lrp6",
    "Ccl12^Ccr2"
]



sasp_top_pairs = [
    "Ceacam1^Havcr2", 
    "Vtcn1^Pdcd1",     ## new 20260120
    "Cd274^Pdcd1",      
    "Cd86^Ctla4",     
    "Cd80^Ctla4",    
    "Cd28^Cd86", 
    "Il1b^Adrb2", 
    "Tnf^Tnfrsf1b",  
    "Cxcl10^Cxcr3",   
    "Ccl4^Ccr5",    
    "Ccl12^Ccr5",
    "Cxcl12^Itga4",  
    "Tnfsf10^Ripk1",   
    "Tgfb1^Cxcr4",
    "Thbs4^Cd47",
    "Thbs4^Cd47"     
]

selected_pairs = [
    "Thbs4^Cd47",    # Cd47: Do not eat me
    "Col1a1^Itga2",
    "Vtcn1^Pdcd1"
]

df = df[df['interaction'].isin(sasp_top_pairs)]

df = df[df['interaction'].isin(region_top_pairs)]

df = df[df['interaction'].isin(selected_pairs)]
df = df[df['cell_pair'].isin(pairs)]
df.to_csv('./results/BR/BR_cluster_cellphone_db_plot.csv')

fig = li.pl.dotplot_by_sample(
    adata=None,
    liana_res=df,
    colour='lr_means',
    size = 'cellphone_pvals',
    inverse_size=True,
    sample_key='cluster',
    figure_size=(8, 25)
)
fig.save('./figures/BR/BR_cellphone_db_dotplot_selected_pairs.png', dpi=300, bbox_inches='tight')


fig = li.pl.dotplot_by_sample(
    adata=None,
    liana_res=df,
    colour='lr_means',
    size = 'cellphone_pvals',
    inverse_size=True,
    sample_key='cluster',
    figure_size=(8, 25)
)
fig.save('./figures/BR/BR_cluster_cellphone_db_dotplot_addVtcn1.png', dpi=300, bbox_inches='tight')

fig = li.pl.dotplot_by_sample(
    adata=None,
    liana_res=df,
    colour='lr_means',
    size = 'cellphone_pvals',
    inverse_size=True,
    sample_key='cluster',
    figure_size=(8, 25)
)
fig.save('./figures/BR/BR_cellphone_db_dotplot_region.png', dpi=300, bbox_inches='tight')
fig.save('./figures/BR/BR_cellphone_db_dotplot_region.pdf', dpi=300, bbox_inches='tight')



fig.save('./figures/BR/BR_cluster_cellphone_db_dotplot——large_b_dc.pdf', dpi=300, bbox_inches='tight')
fig.save('./figures/BR/BR_cluster_cellphone_db_dotplot.svg', dpi=300, bbox_inches='tight')

for region in regions:
    sub_data = df[df['cluster'] == region]
    fig = li.pl.dotplot(adata = None,
              liana_res = sub_data,
              colour='lr_means',
              size='cellphone_pvals',
              inverse_size=True, 
              orderby = 'cellphone_pvals',
              orderby_ascending = True,
              figure_size=(10, 5)
             )
    fig.save('./figures/BR/BR_cluster_'+region+'_cellphone_db_dotplot.pdf', dpi=300, bbox_inches='tight')
    fig.save('./figures/BR/BR_cluster_'+region+'_cellphone_db_dotplot.png', dpi=300, bbox_inches='tight')


