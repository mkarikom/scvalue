import numpy as np
import pandas as pd
import scanpy as sc
from time import time
import matplotlib.pyplot as plt
import os

from scsampler import scsampler
from geosketch import gs
from baselines.treehopper.hoppers import hopper, treehopper, PCATreePartition
from baselines.kernel_herding.scripts.core.model import *

from scvalue import SCValue

from scipy.spatial.distance import directed_hausdorff
from gini import compute_gini
from scipy.spatial import distance_matrix

def make_folders(dataset):
    os.makedirs('experiments/%s' % (dataset), exist_ok=True)

def compute_metrics(adata, adata_sub, cell_type_key='cell_type'):
    original_matrix = adata.obsm['X_pca']
    sketched_matrix = adata_sub.obsm['X_pca']
    hausdorff_distance = directed_hausdorff(original_matrix, sketched_matrix)[0]
    gini_coef = compute_gini(adata, adata_sub, cell_type_key)
    return hausdorff_distance, gini_coef

def write_metric_df(dataset, sketch_method, sketch_size_list, time_list, hausdorff_list, gini_list):
    metric_df = pd.DataFrame({
        'sketch_size': sketch_size_list,
        'sketch_time': time_list,
        'hausdorff_distance': hausdorff_list,
        'gini_coefficient': gini_list
    })
    metric_df.to_csv('experiments/%s/%s.metrics.csv' % (dataset, sketch_method), header=True, index=False)


def generate_outputs(adata, sketch_size_list, sketch_index_list, time_list, dataset, sketch_method, cell_type_key):
    hausdorff_list = []
    gini_list = []
    for i in range(len(sketch_index_list)):
        sketch_index = sketch_index_list[i]
        sketch_size = sketch_size_list[i]
        adata_sub = adata[adata.obs_names[sketch_index]].copy()
        adata_sub.write_h5ad('experiments/%s/%s.%d.h5ad' % (dataset, sketch_method, sketch_size))
        print(adata_sub)

        hausdorff_distance, gini_coef = compute_metrics(adata, adata_sub, cell_type_key)
        hausdorff_list.append(hausdorff_distance)
        gini_list.append(gini_coef)

    write_metric_df(dataset, sketch_method, sketch_size_list, time_list, hausdorff_list, gini_list)

def generate_outputs2(adata, sketch_size_list, time_list, dataset, sketch_method, cell_type_key):
    hausdorff_list = []
    gini_list = []
    for sketch_size in sketch_size_list: 
        adata_sub = sc.read_h5ad('experiments/%s/%s.%d.h5ad' % (dataset, sketch_method, sketch_size))
        print(adata_sub)
        
        hausdorff_distance, gini_coef = compute_metrics(adata, adata_sub, cell_type_key)
        hausdorff_list.append(hausdorff_distance)
        gini_list.append(gini_coef)

    write_metric_df(dataset, sketch_method, sketch_size_list, time_list, hausdorff_list, gini_list)

def eval_geosketch(adata, sketch_size_list, dataset, cell_type_key='cell_type'):
    pca_matrix = adata.obsm['X_pca']
    time_list = []
    sketch_index_list = []
    for sketch_size in sketch_size_list: 
        start = time()
        sketch_index = gs(pca_matrix, sketch_size, replace=False)
        end = time()
        
        sketch_time = end - start
        time_list.append(sketch_time)

        sketch_index_list.append(sketch_index)

    generate_outputs(adata, sketch_size_list, sketch_index_list, time_list, dataset, sketch_method='GeoSketch', cell_type_key=cell_type_key)

def eval_hopper(adata, sketch_size_list, dataset, cell_type_key='cell_type', max_partition_size=1000): # tree hopper actually due to high computational complexity of hopper
    pca_matrix = adata.obsm['X_pca']
    time_list = []
    sketch_index_list = []
    for sketch_size in sketch_size_list: 
        start = time()
        th = treehopper(pca_matrix, partition=PCATreePartition, max_partition_size=max_partition_size)
        th.hop(sketch_size)
        sketch_index = th.path[:sketch_size]
        end = time()
        
        sketch_time = end - start
        time_list.append(sketch_time)

        sketch_index_list.append(sketch_index)

    generate_outputs(adata, sketch_size_list, sketch_index_list, time_list, dataset, sketch_method='Hopper', cell_type_key=cell_type_key)

def get_gamma_range(X): # for KH
    inds = np.random.choice(X.shape[0], size=1000)
    distances = distance_matrix(X[inds], X[inds])
    gamma_0 = np.median(distances)
    gammas = [gamma_0 / i for i in (4, 3, 2, 1, 0.5, 0.33, 0.25, 0.1)]
    return gammas

def eval_kh(adata, sketch_size_list, dataset, cell_type_key='cell_type'):
    pca_matrix = adata.obsm['X_pca']
    time_list = []
    sketch_index_list = []
    for sketch_size in sketch_size_list:
        start = time()
        gammas = get_gamma_range(pca_matrix)
        scale_factor = 1.0 # default
        gamma = gammas[3] * scale_factor
        phi = random_feats(pca_matrix, gamma)
        sketch_index, _, _ = kernel_herding_main(pca_matrix, phi, sketch_size)
        end = time()
        
        sketch_time = end - start
        time_list.append(sketch_time)

        sketch_index_list.append(sketch_index)

    generate_outputs(adata, sketch_size_list, sketch_index_list, time_list, dataset, sketch_method='KH', cell_type_key=cell_type_key)
    
def eval_sphetcher(adata, sketch_size_list, dataset, cell_type_key='cell_type'):
    sketch_index_list = []
    time_list = []
    for sketch_size in sketch_size_list: 
        sketch_indicator_df = pd.read_csv('experiments/%s/Sphetcher_indicator.%d.csv' % (dataset, sketch_size), header=None)
        sketch_index = sketch_indicator_df[sketch_indicator_df[0]==1].index
        sketch_index_list.append(sketch_index)

        time_df = pd.read_csv('experiments/%s/Sphetcher.%d.time.csv' % (dataset, sketch_size), header=0)
        time_list.append(time_df['sketch_time'].values[0])

    generate_outputs(adata, sketch_size_list, sketch_index_list, time_list, dataset, sketch_method='Sphetcher', cell_type_key=cell_type_key)

def eval_scsampler(adata, sketch_size_list, dataset, random_split=None, cell_type_key='cell_type'):
    sketch_method = 'scSampler'
    time_list = []
    for sketch_size in sketch_size_list: 
        start = time()
        adata_sub = scsampler(adata, n_obs=sketch_size, copy=True, random_split=random_split)
        end = time()
        sketch_time = end - start
        time_list.append(sketch_time)

        adata_sub.write_h5ad('experiments/%s/%s.%d.h5ad' % (dataset, sketch_method, sketch_size))

    generate_outputs2(adata, sketch_size_list, time_list, dataset, sketch_method=sketch_method, cell_type_key=cell_type_key)   

def eval_uniform(adata, sketch_size_list, dataset, cell_type_key='cell_type', seed=42):
    time_list = []
    for sketch_size in sketch_size_list: 
        start = time()
        adata_sub = adata.copy()
        sc.pp.subsample(adata_sub, n_obs=sketch_size, random_state=seed)
        end = time()
        sketch_time = end - start
    
        time_list.append(sketch_time)
        adata_sub.write_h5ad('experiments/%s/Uniform.%d.h5ad' % (dataset, sketch_size))

    generate_outputs2(adata, sketch_size_list, time_list, dataset, sketch_method='Uniform', cell_type_key=cell_type_key)   
    
def eval_scvalue(adata, sketch_size_list, dataset, cell_type_key='cell_type', n_trees=100, prop_sampling=False, write_dv=False, strategy='FB', use_rep='X_pca'):
    sketch_method = 'scValue'
    time_list = []
    for sketch_size in sketch_size_list:
        start = time()
        scv = SCValue(adata=adata, sketch_size=sketch_size, cell_type_key=cell_type_key, n_trees=n_trees, prop_sampling=prop_sampling, write_dv=write_dv, strategy=strategy, use_rep=use_rep)
        scv.value()
        adata_sub = scv.adata_sub
        end = time()
        sketch_time = end - start
    
        time_list.append(sketch_time)
        adata_sub.write_h5ad('experiments/%s/%s.%d.h5ad' % (dataset, sketch_method, sketch_size))

    generate_outputs2(adata, sketch_size_list, time_list, dataset, sketch_method=sketch_method, cell_type_key=cell_type_key)   
    
def plot_umap(sketch_size_list, dataset, sketch_method, legend_loc=None, cell_type_key='cell_type', point_size=5):
    for sketch_size in sketch_size_list:
        adata = sc.read_h5ad('experiments/%s/%s.%d.h5ad' % (dataset, sketch_method, sketch_size))

        sc.pl.umap(adata, 
           color=cell_type_key,
           legend_loc=legend_loc,
           show=False,
           size=point_size,
           add_outline=True,
           outline_color=('gray', 'white'),
           frameon=True,
           palette=['#9CAECE', '#ECD7AB', '#E79DB6', '#9EB69D'],
           title='',)
    
        plt.gcf().set_size_inches(4, 4)
        ax = plt.gca()

        color = 'dimgray'
        ax.spines['left'].set_color(color)
        ax.spines['right'].set_color(color)
        ax.spines['top'].set_color(color)
        ax.spines['bottom'].set_color(color)

        ax.xaxis.label.set_color(color)
        ax.yaxis.label.set_color(color)

        ax.tick_params(axis='x', colors=color)
        ax.tick_params(axis='y', colors=color)

        os.makedirs('figures/%s' % (dataset), exist_ok=True)
        plt.savefig('figures/%s/%s.%d.png' % (dataset, sketch_method, sketch_size), dpi=300)
