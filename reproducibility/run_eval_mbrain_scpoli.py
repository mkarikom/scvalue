import numpy as np
import scanpy as sc
import pandas as pd
from scarches.models.scpoli import scPoli
import os
import torch
import pytorch_lightning as pl 

import warnings
warnings.filterwarnings('ignore')

orig_source_adata = sc.read('data/ref.h5ad')
target_adata = sc.read('data/query.h5ad')
target_adata = target_adata[:, orig_source_adata.var_names]

early_stopping_kwargs = {
    "early_stopping_metric": "val_prototype_loss",
    "mode": "min",
    "threshold": 0,
    "patience": 5,
    "reduce_lr": True,
    "lr_patience": 3,
    "lr_factor": 0.1,
}

condition_key = 'study'
cell_type_key = 'cell_type'

path = 'experiments/mbrain/'

def eval_ctp(methods, pretrain_epoch, finetune_epoch, sketch_size, seed, target_adata=target_adata, path=path):

    np.random.seed(seed)
    torch.manual_seed(seed) 
    if torch.cuda.is_available(): 
        torch.cuda.manual_seed_all(seed)
    pl.seed_everything(seed)
    
    n_epochs = pretrain_epoch + finetune_epoch

    acc_list = []
    for method in methods:
        print('scPoli on', method)

        source_adata = sc.read_h5ad(path + method + '.%d.h5ad' % (sketch_size))
        source_adata.X = source_adata.X.toarray()
        source_adata.X = source_adata.X.astype(np.float32)

        scpoli_model = scPoli(
            adata=source_adata,
            condition_keys=condition_key,
            cell_type_keys=cell_type_key,
            embedding_dims=5,
            recon_loss='nb',
        )

        scpoli_model.train(
            n_epochs=n_epochs,
            pretraining_epochs=pretrain_epoch,
            early_stopping_kwargs=early_stopping_kwargs,
            eta=5,
            lr=0.01
        )

        scpoli_query = scPoli.load_query_data(
            adata=target_adata,
            reference_model=scpoli_model,
            labeled_indices=[],
        )

        scpoli_query.train(
            n_epochs=10,#n_epochs,
            pretraining_epochs=8,#pretrain_epoch,
            eta=10,
            lr=0.01
        )

        target_adata_temp = target_adata.copy()
        target_adata_temp.X = target_adata_temp.X.toarray()
        target_adata_temp.X = target_adata_temp.X.astype(np.float32)
        results_dict = scpoli_query.classify(target_adata_temp, scale_uncertainties=True)

        y_true=target_adata_temp.obs[cell_type_key]
        y_pred=results_dict[cell_type_key]["preds"]
        acc = np.mean(y_true == y_pred)
        print("Acc: {}".format(acc))
        acc_list.append(acc)
    
    df = pd.DataFrame({'method': methods, 'acc': acc_list})

    folder = 'experiments/seed%d' % (seed)
    if not os.path.exists(folder):
        os.makedirs(folder)

    df.to_csv('experiments/seed%d/mbrain_acc.%d.csv' % (seed, sketch_size), header=True, index=False)

sketch_percent_list = [0.02, 0.04, 0.06, 0.08, 0.1]
sketch_size_list = [int(41896 * x) for x in sketch_percent_list]

methods = ['scValue', 'Uniform', 'KH', 'GeoSketch', 'Sphetcher', 'Hopper',  'scSampler'] 

pretrain_epoch_list = [4, 4, 4, 8, 8]
finetune_epoch_list = [1, 1, 1, 2, 2]

for j in range(10):
    seed = 42 + j
    for i in range(5):
        sketch_size = sketch_size_list[i]
        pretrain_epoch = pretrain_epoch_list[i]
        finetune_epoch = finetune_epoch_list[i]
        eval_ctp(methods, pretrain_epoch, finetune_epoch, sketch_size, seed)
