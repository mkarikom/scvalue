# scValue

## Desciption
`scValue` is a Python package designed for efficient value-based subsampling of large scRNA-seq datasets. It prioritises cells of higher value (indicative of greater utility for cell type identification) over cells of lower value and allocates more representation to cell types with greater value variability in the resulting subsamples, which support effective downstream machine/deep learning tasks. 

The package has been implemented using `Python` 3.10, `scikit-learn` 1.5.2, `SciPy` 1.14.0, `NumPy` 1.26.4, `pandas` 1.5.3, and `joblib` 1.4.2. Experiments in the paper were conducted using `Scanpy` 1.10.2 and `AnnData` 0.10.8.

## Installation
The package can be installed using pip:
``` python
pip install scvalue
```

## Quick start
Load necessary modules
```python
import pandas as pd
import numpy as np
import scanpy as sc
from scvalue import SCValue
```
### Read and preprocess input scRNA-seq data 
The example dataset of 89,429 mouse T cells (mTC) can be downloaded from https://datasets.cellxgene.cziscience.com/6c253303-8255-4004-8464-5e90dd5846cb.h5ad in the h5ad format. It has already been log-normalised. Here, we perform principal component analysis (PCA) on the top 3,000 highly variable genes (HVGs), as scValue accepts principal components (PCs) as the default features (`n_comps=50`).
```{python}
adata = sc.read_h5ad('mTC.h5ad')
sc.pp.highly_variable_genes(adata, n_top_genes=3000)
sc.pp.pca(adata, n_comps=50, use_highly_variable=True)
```

### Subsample by scValue
Subsample 10% cells of mTC and return a new AnnData containing the subsample. Alternatively, `sketch_size` can be specified as an integer, representing the exact number of cells to be subsampled. Users can also specify:
- the type of features `use_rep` (default "X_pca")
- `cell_type_key` (default "cell_type") in the original anndata
- the number of trees, `n_trees`, in the random forest (default 100)
- the subsampling `strategy`, one of the three options, listed in an increasing order of cell type separation: 
    - full binning, "FB" (default)
    - mean-threshold binning, "MTB"
    - top-pick, "TP"

Instead of value-weighted subsampling, proportional subsampling, `prop_sampling`, can be carried out if the original cell types are highly balanced or if users wish to maintain the original cell type proportions. In addition, the data value of individual cells can be saved to disk for future reference by setting `write_dv=True`. Lastly, a random `seed` (default 42) can be specified for random forest fitting. 
```{python}
# create an instance of the SCValue class
scv = SCValue(adata=adata, 
              sketch_size=0.1, 
              use_rep="X_pca",
              cell_type_key="cell_type", 
              n_trees=100,
              strategy="FB", 
              prop_sampling=False,
              write_dv=False,
              seed=42)
scv.value() # perform value computation and value-based subsampling
adata_sub = scv.adata_sub
```

`adata_sub` can be plotted using `plot_umap` in `reproducibility/eval_time_hd_gini.py` to reproduce the demo (Section 3.1) in the paper. All datasets used in the paper are publicly available, as detailed in Supplementary Section S1.

## Reference
Li Huang & Dongsheng Chen, scValue: value-based subsampling of large-scale single-cell transcriptomic data. 2024.
