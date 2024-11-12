import numpy as np
import pandas as pd
import scanpy as sc

from eval_time_hd_gini import *

# Take scArches' pbmc data as an example
adata = sc.read_h5ad('scarches_data/pbmc_source.h5ad')
cell_type_column = 'CellType'
print(adata)
print(adata.obs[cell_type_column].unique())

dataset = 'scarches_pbmc'

sketch_percent_list = [0.02, 0.04, 0.06, 0.08, 0.1]
sketch_size_list = [int(adata.shape[0] * x) for x in sketch_percent_list]
    
eval_scvalue(adata, sketch_size_list, dataset=dataset, cell_type_column=cell_type_column)

eval_scsampler(adata, sketch_size_list, dataset, random_split=4, cell_type_column=cell_type_column)

eval_kh(adata, sketch_size_list, dataset, cell_type_column=cell_type_column)

eval_hopper(adata, sketch_size_list, dataset, cell_type_column=cell_type_column)

eval_geosketch(adata, sketch_size_list, dataset, cell_type_column=cell_type_column)

eval_uniform(adata, sketch_size_list, dataset, cell_type_column=cell_type_column)

eval_sphetcher(adata, sketch_size_list, dataset, cell_type_column=cell_type_column)
