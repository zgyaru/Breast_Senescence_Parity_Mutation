import os
import anndata as ad
import matplotlib.pyplot as plt
import scanpy as sc
import scvi
import torch

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
os.environ["SLURM_NTASKS_PER_NODE"] = "16"


def run_scVI(adata, batchID='batch', n_dims=20, n_layers=2, max_epochs=400):
    train_adata = adata[:, adata.var['highly_variable']].copy()
    scvi.model.SCVI.setup_anndata(train_adata, batch_key=batchID)

    arches_params = {
        "use_layer_norm": "both",
        "use_batch_norm": "none",
        "encode_covariates": True,
        "dropout_rate": 0.2,
        "n_layers": n_layers,
    }

    vae = scvi.model.SCVI(
        train_adata,
        n_latent=n_dims,
        gene_likelihood="nb",
        **arches_params
    )

    vae.train(
        early_stopping=True,
        train_size=0.9,
        early_stopping_patience=45,
        max_epochs=max_epochs,
        batch_size=1024,
        accelerator="cpu",
        limit_train_batches=20
    )

    adata.obsm['X_scVI'] = vae.get_latent_representation()
    return adata, vae


def post_dimred_processing(adata, dimred='X_scVI', run_clustering=True):
    sc.pp.neighbors(adata, use_rep=dimred, n_neighbors=15)
    sc.tl.umap(adata)

    if run_clustering:
        for res in [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
            sc.tl.leiden(adata, resolution=res, key_added=f"leiden_{str(res).replace('.', '_')}")

    sc.tl.diffmap(adata)
    return adata


def prepare_subset(adata, compartment, hvg_n=5000):
    sub = adata[adata.obs['ct_level0'] == compartment].copy()
    sub.varm = None
    sub.obsp = None
    sub.obsm = None
    sub.X = sub.layers['log1p_norm'].copy()
    sc.pp.highly_variable_genes(sub, n_top_genes=hvg_n)
    sub.X = sub.layers['raw'].copy()
    return sub


save_dir = '../integration/scvi_output/'
adata = ad.read_h5ad('../integration/Level1/adata_full_scvi_level1.h5ad')

os.makedirs(save_dir + 'loss', exist_ok=True)
os.makedirs(save_dir + 'umap', exist_ok=True)
os.makedirs(save_dir + 'adata', exist_ok=True)
os.makedirs(save_dir + 'model', exist_ok=True)

compartments = {
    'Epithelial': 'epi',
    'Lymphoid': 'lym',
    'Myeloid': 'mye',
    'Stroma': 'stromal'
}

for compartment, prefix in compartments.items():
    print(f'Processing {compartment}...')

    sub = prepare_subset(adata, compartment, hvg_n=5000)
    sub, vae = run_scVI(sub, batchID='Batch', n_dims=30, n_layers=2, max_epochs=400)

    plt.clf()
    plt.plot(vae.history["elbo_train"], label="train")
    plt.plot(vae.history["elbo_validation"], label="validation")
    plt.legend()
    plt.savefig(save_dir + f'loss/elbo_{prefix}.png')

    plt.clf()
    plt.plot(vae.history["reconstruction_loss_train"], label="train")
    plt.plot(vae.history["reconstruction_loss_validation"], label="validation")
    plt.legend()
    plt.savefig(save_dir + f'loss/reconstruction_{prefix}.png')

    sc.settings.figdir = save_dir
    sc.pp.neighbors(sub, use_rep='X_scVI', n_neighbors=15)
    sc.tl.umap(sub)
    sc.pl.umap(sub, color='Batch', save=f'/{prefix}_batch.png')

    sub = post_dimred_processing(sub, dimred='X_scVI', run_clustering=True)
    sub.write_h5ad(save_dir + f'adata/{prefix}_scvi.h5ad')
    vae.save(save_dir + f'model/{prefix}_model', overwrite=True)