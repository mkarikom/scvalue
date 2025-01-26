import scanpy as sc
import celltypist
import numpy as np
import pandas as pd
import os

adata_test = sc.read_h5ad('data/test.h5ad')
path = 'experiments/cxg_min/'

def eval_ctp(methods, sketch_size, seed, adata_test=adata_test, path=path):

    np.random.seed(seed)

    acc_list = []
    for method in methods:
        print('celltypist on', method)
        adata_train = sc.read_h5ad(path + method + '.%d.h5ad' % (sketch_size))
        print(adata_train)

        model_fs = celltypist.train(adata_train, 'cell_type', n_jobs = 10, max_iter = 5, use_SGD = True)
        gene_index = np.argpartition(np.abs(model_fs.classifier.coef_), -100, axis = 1)[:, -100:]
        gene_index = np.unique(gene_index)
        print(f"Number of genes selected: {len(gene_index)}")

        model = celltypist.train(adata_train[:, gene_index], 'cell_type', check_expression = False, n_jobs = 10, max_iter = 10, use_SGD = True)
        model.write('data/model_train_sketch.pkl')

        predictions = celltypist.annotate(adata_test, model = 'data/model_train_sketch.pkl')
        print(predictions.predicted_labels)

        acc = np.mean(adata_test.obs['cell_type'].astype(str) == predictions.predicted_labels['predicted_labels'].astype(str))
        print("Test Acc: {}".format(acc))
        acc_list.append(acc)
    
    df = pd.DataFrame({'method': methods, 'acc': acc_list})

    folder = 'experiments/seed%d' % (seed)
    if not os.path.exists(folder):
        os.makedirs(folder)

    df.to_csv('experiments/seed%d/cxg_min_acc.%d.csv' % (seed, sketch_size), header=True, index=False)


sketch_percent_list = [0.02, 0.04, 0.06, 0.08, 0.1]
sketch_size_list = [int(32768 * x) for x in sketch_percent_list]

methods = ['scValue', 'Uniform', 'GeoSketch', 'Sphetcher', 'Hopper', 'KH', 'scSampler'] 

for j in range(10):
    seed = 42 + j
    for i in range(5):
        sketch_size = sketch_size_list[i]
        eval_ctp(methods, sketch_size, seed)
