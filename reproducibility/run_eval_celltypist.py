import numpy as np
import scanpy as sc
import celltypist
import time

np.random.seed(42)

adata_James = sc.read_h5ad('celltypist_data/Colon_cell_atlas.h5ad')
print(adata_James.obs.cell_type.unique())

orders = [
    'Activated CD4 T',
    'Th1',
    'Tfh',
    'CD8 T',
    'cycling gd T',
    'Tcm',
    'gd T',
    'Th17',
    'Treg'
]

path = 'experiments/celltypist_gut/'
methods = ['Uniform', 'GeoSketch', 'Sphetcher', 'Hopper', 'KH', 'scSampler', 'scValue']
sketch_size = 5485 # 10% of reference data

for method in methods:
    print('Celltypist on', method)
    adata = sc.read_h5ad(path + '%s.%d.h5ad' % (method, sketch_size))
    print(adata)

    t_start = time.time()
    model = celltypist.train(adata, 'Integrated_05', check_expression = False, n_jobs = 10, max_iter = 100)
    t_end = time.time()
    print(f"Time elapsed: {t_end - t_start} second")

    # Save the model.
    model.write('celltypist_data/model_from_Elmentaite_2021.%s.pkl' % (method))

    t_start = time.time()
    predictions = celltypist.annotate(adata_James, model = 'celltypist_data/model_from_Elmentaite_2021.%s.pkl' % (method), majority_voting = True)
    t_end = time.time()
    print(f"Time elapsed: {t_end - t_start} seconds")

    celltypist.dotplot(predictions, use_as_reference = 'cell_type', use_as_prediction = 'majority_voting', 
                       show=False, reference_order=orders, save='%s_T.png' % (method), title='')

    adata_pred = predictions.to_adata()
    adata_pred.write_h5ad(path + method + '.pred.h5ad')
