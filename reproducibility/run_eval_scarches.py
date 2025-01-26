import scanpy as sc
import torch
import scarches as sca
from scarches.dataset.trvae.data_handling import remove_sparsity
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=UserWarning)

torch.set_printoptions(precision=3, sci_mode=False, edgeitems=7)

target_adata = sc.read_h5ad('scarches_data/pbmc_target.h5ad')

path = 'experiments/scarches_pbmc/'
methods = ['Uniform', 'GeoSketch', 'Sphetcher', 'Hopper', 'KH', 'scSampler', 'scValue'] 

condition_key = 'Method'
cell_type_key = 'CellType'

def eval_scarches(methods, sketch_size, seed):
    np.random.seed(seed)
    torch.manual_seed(seed)

    acc_list = []
    for method in methods:
        print('scarches on', method)
        source_adata = sc.read_h5ad(path + method + '.%d.h5ad' % (sketch_size))
        print(source_adata)

        sca.models.SCVI.setup_anndata(source_adata, batch_key=condition_key, labels_key=cell_type_key)

        vae = sca.models.SCVI(
            source_adata,
            n_layers=2,
            encode_covariates=True,
            deeply_inject_covariates=False,
            use_layer_norm="both",
            use_batch_norm="none",
        )

        vae.train(max_epochs=100)
        
        scanvae = sca.models.SCANVI.from_scvi_model(vae, unlabeled_category = "Unknown")

        print("Labelled Indices: ", len(scanvae._labeled_indices))
        print("Unlabelled Indices: ", len(scanvae._unlabeled_indices))

        scanvae.train(max_epochs=10)

        reference_latent = sc.AnnData(scanvae.get_latent_representation())
        reference_latent.obs["cell_type"] = source_adata.obs[cell_type_key].tolist()
        reference_latent.obs["batch"] = source_adata.obs[condition_key].tolist()

        reference_latent.obs['predictions'] = scanvae.predict()
        print("Train Acc: {}".format(np.mean(reference_latent.obs.predictions == reference_latent.obs.cell_type)))

        ref_path = 'scarches_data/'
        scanvae.save(ref_path, overwrite=True)

        model = sca.models.SCANVI.load_query_data(
            target_adata,
            ref_path,
            freeze_dropout = True
        )

        model._unlabeled_indices = np.arange(target_adata.n_obs)
        model._labeled_indices = []
        print("Labelled Indices: ", len(model._labeled_indices))
        print("Unlabelled Indices: ", len(model._unlabeled_indices))

        model.train(
            max_epochs=10,
            plan_kwargs=dict(weight_decay=0.0),
            check_val_every_n_epoch=10,
        )

        query_latent = sc.AnnData(model.get_latent_representation())
        query_latent.obs['cell_type'] = target_adata.obs[cell_type_key].tolist()
        query_latent.obs['batch'] = target_adata.obs[condition_key].tolist()

        #surgery_path = 'surgery_model'
        #model.save(surgery_path, overwrite=True)

        query_latent.obs['predictions'] = model.predict()
        acc = np.mean(query_latent.obs.predictions == query_latent.obs.cell_type)
        print("Test Acc: {}".format(acc))
        acc_list.append(acc)

    df = pd.DataFrame({'method': methods, 'acc': acc_list})

    folder = 'experiments/seed%d' % (seed)
    if not os.path.exists(folder):
        os.makedirs(folder)
    df.to_csv('experiments/seed%d/acc.%d.csv' % (seed, sketch_size), header=True, index=False)

sketch_percent_list = [0.02, 0.04, 0.06, 0.08, 0.1]
sketch_size_list = [int(16941 * x) for x in sketch_percent_list]
for i in range(10):
    seed = 42 + i
    for sketch_size in sketch_size_list:
        eval_scarches(methods, sketch_size, seed)
