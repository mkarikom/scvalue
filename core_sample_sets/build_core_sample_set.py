import scanpy as sc
import numpy as np
import pandas as pd
import celltypist
import os

def pred_by_celltypist(adata, # the scRNA-seq dataset (anndata format) used to build the core sample set
                       pred_folder, # the output folder where CellTypist predictions are saved
                       model_path, # the path of the CellTypist model used to predict cell types of the dataset and compute the corresponding confidence scores
                       cell_type_key='cell_type'
                       ):

    """
    Compute confidence scores ([0, 1]) of cell types in a given scRNA-seq dataset (anndata format) by either: 
     (1) an established CellTypist model obtained from the official model repository https://www.celltypist.org/models;
     (2) a custom CellTpyist model trained with high-quality data by the user.
    The model should be in the pickle format, such as "Immune_All_High.pkl".
    """

    os.makedirs(pred_folder, exist_ok=True)

    # check if the expression matrix is log1p normalized to 10000 counts per cell, as required by CellTypist to give reliable predictions
    print('log1p normalized to 10000 counts per cell: ', np.expm1(adata.X).sum(axis=1))
    predictions = celltypist.annotate(adata, model = model_path, majority_voting = True)
    
    adata_pred = predictions.to_adata()
    # the metadata of all cells, together with their predicted cell types and confidence scores
    adata_pred.obs.to_csv(pred_folder + '/obs_pred.csv', header=True, index=True)

    pred_df = adata_pred.obs[[cell_type_key, 'majority_voting']].copy()
    pred_df.drop_duplicates(inplace=True)
    # Unique combinations of original cell types and predicted cell types
    pred_df.to_csv(pred_folder + '/pred_summary.csv', header=True, index=False)

    return pred_df

def build_core_sample_set(adata, # the scRNA-seq dataset (anndata format) used to build the core sample set
                          pred_kept_path, # the user needs to specify which combinations of original and predicted cell types are matched in the pred_summary.csv,
                                          # delete unmatched ones, and keep only these matched as input to build the core sample set
                          pred_folder, # the pred folder where the metadata, predictions, and confidence scores of all cells have been saved                           
                          data_folder, # the output folder where the core sample set (anndata) is saved
                          cell_type_key='cell_type', 
                          conf_score_threshold=0.5 # cells with matched original and predicted cell types and confidence scores greater than the theshold are included in the core sample set
                          ):

    """
    In the given scRNA-seq dataset (anndata format), cells that meet both criteria are included in the core sample set: 
     (1) original and predicted cell types are matched, either exactly the same type or belonging to the same broader category 
        (unassigned/unknown cells can also be included);
     (2) confidence scores are above the specified theshold, default 0.5.
    """

    os.makedirs(data_folder, exist_ok=True)

    obs_pred = pd.read_csv(pred_folder + '/obs_pred.csv', header=0, index_col = 0)
    pred_kept = pd.read_csv(pred_kept_path, header=0)

    pred_kept['key'] = list(zip(pred_kept[cell_type_key], pred_kept['majority_voting']))
    obs_pred['key'] = list(zip(obs_pred[cell_type_key], obs_pred['majority_voting']))

    obs_match = obs_pred[obs_pred['key'].isin(pred_kept['key'])].drop(columns='key')
    obs_core = obs_match[obs_match['conf_score'] >= conf_score_threshold].copy()

    adata_core = adata[obs_core.index.astype(str), :].copy()
    adata_core.write_h5ad(data_folder + '/pbmc_source_core.h5ad')
    print(adata.obs[cell_type_key].value_counts())
    print(adata_core.obs[cell_type_key].value_counts())

    return adata_core
