import numpy as np
import pandas as pd

def compute_gini(adata, adata_sub, cell_type_column):
    ori_prop = adata.obs[cell_type_column].value_counts()/adata.n_obs
    ori_prop = pd.DataFrame(ori_prop)
    ori_prop.columns = ['count']

    sub_prop = adata_sub.obs[cell_type_column].value_counts()/adata_sub.n_obs
    sub_prop = pd.DataFrame(sub_prop)
    sub_prop.columns = ['count']
    
    prop_df = pd.merge(ori_prop, sub_prop, how='left', left_index = True, right_index = True)
    prop_df = prop_df.fillna(0)

    print(prop_df)

    gini_res = gini(prop_df['count_y'].values)

    return(gini_res)

def gini(array):
    """Calculate the Gini coefficient of a numpy array."""
    # based on bottom eq:
    # http://www.statsdirect.com/help/generatedimages/equations/equation154.svg
    # from:
    # http://www.statsdirect.com/help/default.htm#nonparametric_methods/gini.htm
    # All values are treated equally, arrays must be 1d:
    array = array.flatten()
    if np.amin(array) < 0:
        # Values cannot be negative:
        array -= np.amin(array)
    # Values cannot be 0:
    array = array + 0.0000001
    # Values must be sorted:
    array = np.sort(array)
    # Index per array element:
    index = np.arange(1,array.shape[0]+1)
    # Number of array elements:
    n = array.shape[0]
    # Gini coefficient:
    return ((np.sum((2 * index - n  - 1) * array)) / (n * np.sum(array)))