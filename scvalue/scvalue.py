import numpy as np
import pandas as pd
from numbers import Integral, Real
from typing import Optional, Union

import scipy
from scipy.sparse import issparse

from sklearn.exceptions import DataConversionWarning
from sklearn.tree._tree import DTYPE, DOUBLE
from sklearn.ensemble import RandomForestClassifier
from sklearn.utils import check_random_state, compute_sample_weight
from sklearn.utils.parallel import delayed
from sklearn.utils.multiclass import type_of_target
from sklearn.utils.validation import _check_sample_weight

from joblib import Parallel

import warnings
from warnings import catch_warnings, simplefilter, warn
warnings.filterwarnings("ignore", category=UserWarning, module="sklearn")

MAX_INT = np.iinfo(np.int32).max


class SCValue:
    """
    Data value-based subsampling of large-scale single-cell transcriptomic (scRNA-seq) data.

    This class is instantiated with an AnnData object, which contains large scRNA-seq data 
    and cell type information. Two steps are taken to obtain a subsample (sketch) of the data:
    Step 1:
        Computation of each cell's data values (DV) in terms of out-of-bag (OOB) estimate by 
        instantiating the nested class RandomForestSCV and fitting a random forest (RF) model. 
        Here, DV of a cell is defined as how helpful (high DV) or harmful (low DV) the cell is
        in learning the differences of cell types in RF.
    Step 2:
        DV-weighted sampling is carried out to determine the number of cells to be subsampled in 
        each cell type, aiming to enhance the cell type balanceness in the sketch. Then, three 
        cell selection strategies are available: full binning (FB), mean-threshold binning (MTB),
        or top-pick (TP). Alternatively, proportional sampling can be conducted to respect the 
        original type distribution.
    
    Paramters
    ----------
    sketch_size : Union[int, float], default=None
        The number or percentage of cells to be subsampled from AnnData.
    
    use_rep : str, default='X_rep'
        The cell representations for fitting the random forest and computing DVs.

    cell_type_key : str, default='cell_type'
        The name of the cell type column in AnnData.obs.
    
    n_trees : int, default=100
        The number of trees in the random forest.

    strategy : {'FB', 'MTB', 'TP'}, default='FB'
        The cell selection strategies after determination of the cell number in each type.
        1)  'FB' : full binning
            For cells of each type, their DVs are divided into bins with 0.1 intervals ranging
            from 0.0 to 1.0. Within each bin, cells of the highest DVs are selected to construct
            the sketch.
        2)  'MTB' : 'mean-threshold binning'
            The same binning as FB is performed; only those bins above the DV average are considered.
        3)  'TP' : top-pick
            No binning is involved. Simply select cells with the highest DVs in each type.

    prop_sampling : bool, default=False
        Whether to perform proportional sampling (respecting the original cell type distribution) 
        or DV-weighted sampling (addressing the cell type imbalanceness in the original data).
    
    write_dv : bool, default=False
        Whether to write the computed DVs for all cells in the original AnnData to disk.
    
    seed : int, default=42
        To initialize a random number generator (RNG) in random forest. Ensure reproducibility.

    Attributes
    ----------
    rf_model : class
        The instantiated RandomForestSCValue class, containing n_trees and using all available CPU cores.

    dv : pandas.core.frame.DataFrame
        The pandas DataFrame holding DVs and cell types for all cells in the original AnnData.

    adata_sub : AnnData
        The sketch of the original AnnData.
    """
    def __init__(
            self,
            adata=None,
            sketch_size: Optional[Union[int, float]] = None,
            use_rep='X_pca', # can be gene features 'X' or any representations stored in '.obsm'
            cell_type_key='cell_type',
            n_trees=100,
            strategy='FB', # FB: full binning, MTB: mean-threshold binning, TP: top-pick
            prop_sampling=False,
            write_dv=False,
            seed=42,
    ):
        if adata is None:
            raise ValueError('Anndata object is required.')
        
        n_obs = adata.n_obs
        if sketch_size is None:
            raise ValueError('Sketch size should be specified.')
        elif isinstance(sketch_size, int):
            sketch_size = min(n_obs, sketch_size)
        elif isinstance(sketch_size, float):
            if 0 < sketch_size <= 1:
                sketch_size = int(n_obs * sketch_size)
            else:
                raise ValueError('Sketch percentage should be between 0 and 1.')
        else:
            raise ValueError('Sketch size should be an integer or a float.')

        if use_rep != 'X' and use_rep not in adata.obsm:
            raise ValueError(f'{use_rep} does not exist in adata.obsm. Please provide a valid representation.')
        if cell_type_key not in adata.obs.columns:
            raise ValueError(f'{cell_type_key} does not exist in adata.obs. Please provide a valid cell type key.')

        self.adata = adata
        self.sketch_size = sketch_size
        self.use_rep = use_rep
        self.cell_type_key = cell_type_key
        self.n_trees = n_trees
        self.strategy = strategy
        self.prop_sampling = prop_sampling
        self.write_dv = write_dv
        self.seed = seed

        np.random.seed(self.seed)
        self.rf_model = self.RFValue(n_estimators=self.n_trees, n_jobs=-1) 
        
        self.dv = None
        self.adata_sub = None

    def hamilton(self, exact_alloc, size):
        int_alloc = np.floor(exact_alloc).astype(int)
        remaining = size - int_alloc.sum()
        remainders = exact_alloc - int_alloc
        indices = remainders.sort_values(ascending=False).index
        for i in range(remaining):
            int_alloc[indices[i]] += 1
        return int_alloc

    def get_prop(self, df, col, size):
        counts_series = df[col].value_counts()
        exact_alloc = counts_series / counts_series.sum() * size
        return self.hamilton(exact_alloc, size)

    def get_weighted_prop(self, df, size):
        exact_alloc = df['weighted_prop'] * size
        return self.hamilton(exact_alloc, size)

    def train_rf(self):
        if self.use_rep=='X' or (self.use_rep == 'X_pca' and 'X_pca' not in self.adata.obsm):
            print('Use counts as representations for computing dv')
            if isinstance(self.adata.X, scipy.sparse.csr_matrix):
                self.adata.X = self.adata.X.toarray()
            rep_matrix = self.adata.X
        else:
            print(f'Use {self.use_rep} as representations for computing dv')
            rep_matrix = self.adata.obsm[self.use_rep]

        cell_type_mapping = {cell_type: idx for idx, cell_type in enumerate(self.adata.obs[self.cell_type_key].unique())}
        self.adata.obs['cell_type_numeric'] = self.adata.obs[self.cell_type_key].map(cell_type_mapping)
        cell_type_labels = self.adata.obs['cell_type_numeric'].values
        self.rf_model.fit(rep_matrix, cell_type_labels)

        dv_df = pd.DataFrame(self.rf_model.compute_rf_dv(rep_matrix, cell_type_labels), columns=['dv'], index=self.adata.obs_names)
        return dv_df

    def get_dv_bins(self, threshold):
        bins = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        if threshold == 0:
            return bins
        else:
            above = [x for x in bins if x > threshold]
            dv_bins = [0] + [threshold] + above
            return dv_bins

    def reallocate(self, df, excess):
        df['sample_count'] = df.apply(lambda row: row['count'] if row['sample_count'] > row['count'] 
                                      else row['sample_count'], axis=1)
        while excess > 0:
            valid_df = df[df['sample_count'] < df['count']].copy()
            if valid_df.empty:
                break

            valid_df['alloc_prop'] = valid_df['weighted_prop'] / valid_df['weighted_prop'].sum()
            valid_df['receivable'] = valid_df['count'] - valid_df['sample_count']
            valid_df['alloc'] = np.floor((valid_df['alloc_prop'] * excess)).astype(int)

            if valid_df['alloc'].max() == 0:
                indices = valid_df['alloc_prop'].sort_values(ascending=False).index
                for i in range(int(excess)):
                    valid_df.loc[indices[i], 'alloc'] = 1
            else:
                valid_df['alloc'] = valid_df[['alloc', 'receivable']].min(axis=1)
            
            df.loc[valid_df.index, 'sample_count'] += valid_df['alloc']
            excess -= valid_df['alloc'].sum()
        
        df['sample_count'] = df['sample_count'].astype('int')
        return df

    def value(self):
        warnings.filterwarnings("default", category=self.SketchWarning)
        print(f'Start: sketch {self.sketch_size}')
        
        dv_df = self.train_rf()
        self.adata.obs['dv'] = dv_df['dv']
        obs_df = self.adata.obs.copy()
        obs_df = obs_df.sort_values(by='dv', ascending=False)
        info_df = obs_df.groupby(self.cell_type_key, observed=False)['dv'].describe()
        self.dv = self.adata.obs[[self.cell_type_key, 'dv']].copy()

        if self.write_dv:
            self.dv.to_csv('dv_values.csv', header=True, index=True)
            info_df.to_csv('dv_summary.csv', header=True, index=True)

        if self.prop_sampling:
            print('Perform proportional sampling for each cell type')
            sample_cell_type_counts = self.get_prop(obs_df, self.cell_type_key, self.sketch_size)
        else:
            print('Perform weighted sampling for each cell type')
            info_df['cv'] = info_df['std'] / info_df['mean']
            info_df['cv_adj'] = info_df['cv'] * np.sqrt(info_df['count'])
            info_df['expected_count'] = info_df['count'] / info_df['cv_adj']
            
            info_df['expected_count'] = info_df.apply(lambda row: row['count'] if (row['mean']== 0 or row['std']==0) else row['expected_count'], axis = 1)
            info_df['weighted_prop'] = info_df['expected_count'] / info_df['expected_count'].sum()
            info_df['sample_count'] = self.get_weighted_prop(info_df, self.sketch_size)
            
            info_df['excess'] = info_df['sample_count'] - info_df['count']
            info_df['excess'] = info_df['excess'].apply(lambda x: max(x, 0))
            
            excess = info_df['excess'].sum()
            if excess > 0:
                warnings.warn('For each cell type, sketch size cannot be greater than the original size. Reallocating...', self.SketchWarning)
                info_df = self.reallocate(info_df, excess)

            sample_cell_type_counts = info_df['sample_count']

        sampled_indices = []

        if self.strategy != 'TP':
            print('Get bin-based stratified samples for each cell type')
            for cell_type in sample_cell_type_counts.index:
                obs_df_sub = obs_df[obs_df[self.cell_type_key]==cell_type].copy()

                if self.strategy == 'FB':
                    dv_bins = self.get_dv_bins(0)
                elif self.strategy == 'MTB':
                    bin_threshold = np.round(info_df.loc[cell_type, 'mean'] * 10) / 10
                    dv_bins = self.get_dv_bins(bin_threshold)
                else:
                    raise ValueError("Please provide a valid strategy: one of FB, MTB, or TP")
            
                print('Sketch:', cell_type)
                print('dv bins:', dv_bins)

                obs_df_sub['dv_interval'] = pd.cut(obs_df_sub['dv'], bins=dv_bins, include_lowest=True)
                undersample_size = sample_cell_type_counts[cell_type]
                sample_interval_counts = self.get_prop(obs_df_sub, 'dv_interval', undersample_size)
                sample_interval_counts = sample_interval_counts[sample_interval_counts != 0]
            
                print('Sketch size:', undersample_size)
                for interval in sample_interval_counts.index:
                    obs_df_sub_interval = obs_df_sub[obs_df_sub['dv_interval']==interval].copy()
                    interval_undersample_size = sample_interval_counts[interval]
                    sampled_indices.extend(obs_df_sub_interval.index[:interval_undersample_size])
        else:
            print('Get top dv samples for each cell type')
            for cell_type in sample_cell_type_counts.index:
                print('Sketch:', cell_type)
                undersample_size = sample_cell_type_counts[cell_type]
                print('Sketch size:', undersample_size)
                obs_df_sub = obs_df[obs_df[self.cell_type_key]==cell_type].copy()
                sampled_indices.extend(obs_df_sub.index[:undersample_size])
        print(f'Done: sketch {self.sketch_size}')

        self.adata_sub = self.adata[sampled_indices].copy()
        return self.adata_sub

    class SketchWarning(UserWarning):
        pass

    class RFValue(RandomForestClassifier):
        """
        Inherits the RandomForestClassifier class in sklearn (Version 1.5.2) under the path
        scikit-learn/sklearn/ensemble/_forest.py. RFValue also incorporates three functions,
        i.e., _get_n_samples_bootstrap, _generate_sample_indices, and _parallel_build_trees 
        from the same path. These functions are utilized to support tree-fitting and data 
        value computation. 
        
        Implements computation of OOB estimate as data value based on Kwon and Zou's ICML 
        paper in 2023.

        Creates two additional attributes
        1) index_counts : ndarray of shape (n_trees, n_samples)
           The number of sample occurances in each boostrap dataset.
        2) rf_dv : ndarray of shape (n_samples,)
           The computed data value for each sample, i.e., the average OOB accuracy when the 
           sample is not selected in the bootstrap dataset. Please refer to Kwon and Zou's 
           paper for more detail.

        Parameters are the same as those of RandomForestClassifier. Within RFValue, 

        References:
            [1] L. Breiman, "Random Forests", Machine Learning, 45(1), 5-32, 2001.
            [2] Y. Kwon and J. Zou, "Data-oob: Out-of-bag estimate as a simple and efficient 
                data value", International Conference on Machine Learning, 18135-18152, 2023.
        """
        def __init__(
            self,
            n_estimators=100,
            *,
            criterion="gini",
            max_depth=None,
            min_samples_split=2,
            min_samples_leaf=1,
            min_weight_fraction_leaf=0.0,
            max_features="sqrt",
            max_leaf_nodes=None,
            min_impurity_decrease=0.0,
            bootstrap=True,
            oob_score=False,
            n_jobs=None,
            random_state=None,
            verbose=0,
            warm_start=False,
            class_weight=None,
            ccp_alpha=0.0,
            max_samples=None,
            monotonic_cst=None,
        ):
            super().__init__(
                n_estimators=n_estimators,
                criterion = criterion,
                max_depth = max_depth,
                min_samples_split = min_samples_split,
                min_samples_leaf = min_samples_leaf,
                min_weight_fraction_leaf = min_weight_fraction_leaf,
                max_features = max_features,
                max_leaf_nodes = max_leaf_nodes,
                min_impurity_decrease = min_impurity_decrease,
                bootstrap=bootstrap,
                oob_score=oob_score,
                n_jobs=n_jobs,
                random_state=random_state,
                verbose=verbose,
                warm_start=warm_start,
                class_weight=class_weight,
                ccp_alpha = ccp_alpha,
                max_samples=max_samples,
                monotonic_cst=monotonic_cst)

        def compute_rf_dv(self, X, y):
            tree_dv = [dict(zip(np.where(self.index_counts[i] == 0)[0], 
                                       (tree.predict(X[self.index_counts[i] == 0]) == y[self.index_counts[i] == 0]).astype(float))) 
                                       for i, tree in enumerate(self.estimators_)]
            self.rf_dv = pd.DataFrame(tree_dv)[np.arange(len(X))].mean(axis=0).to_numpy()
            return self.rf_dv

        def _get_n_samples_bootstrap(self, n_samples, max_samples):
            """
            Get the number of samples in a bootstrap sample.

            Parameters
            ----------
            n_samples : int
                Number of samples in the dataset.
            max_samples : int or float
                The maximum number of samples to draw from the total available:
                    - if float, this indicates a fraction of the total and should be
                    the interval `(0.0, 1.0]`;
                    - if int, this indicates the exact number of samples;
                    - if None, this indicates the total number of samples.

            Returns
            -------
            n_samples_bootstrap : int
                The total number of samples to draw for the bootstrap sample.
            """
            if max_samples is None:
                return n_samples

            if isinstance(max_samples, Integral):
                if max_samples > n_samples:
                    msg = "`max_samples` must be <= n_samples={} but got value {}"
                    raise ValueError(msg.format(n_samples, max_samples))
                return max_samples

            if isinstance(max_samples, Real):
                return max(round(n_samples * max_samples), 1)

        def _generate_sample_indices(self, random_state, n_samples, n_samples_bootstrap):
            """
            Private function used to _parallel_build_trees function."""

            random_instance = check_random_state(random_state)
            sample_indices = random_instance.randint(
                0, n_samples, n_samples_bootstrap, dtype=np.int32
            )

            return sample_indices

        def _parallel_build_trees(
            self,
            tree,
            bootstrap,
            X,
            y,
            sample_weight,
            tree_idx,
            n_trees,
            verbose=0,
            class_weight=None,
            n_samples_bootstrap=None,
            missing_values_in_feature_mask=None,
        ):
            """
            Private function used to fit a single tree in parallel."""
            if verbose > 1:
                print("building tree %d of %d" % (tree_idx + 1, n_trees))

            if bootstrap:
                n_samples = X.shape[0]
                if sample_weight is None:
                    curr_sample_weight = np.ones((n_samples,), dtype=np.float64)
                else:
                    curr_sample_weight = sample_weight.copy()

                indices = self._generate_sample_indices(
                    tree.random_state, n_samples, n_samples_bootstrap
                )
                sample_counts = np.bincount(indices, minlength=n_samples)
                curr_sample_weight *= sample_counts

                if class_weight == "subsample":
                    with catch_warnings():
                        simplefilter("ignore", DeprecationWarning)
                        curr_sample_weight *= compute_sample_weight("auto", y, indices=indices)
                elif class_weight == "balanced_subsample":
                    curr_sample_weight *= compute_sample_weight("balanced", y, indices=indices)

                tree._fit(
                    X,
                    y,
                    sample_weight=curr_sample_weight,
                    check_input=False,
                    missing_values_in_feature_mask=missing_values_in_feature_mask,
                )
                return tree, curr_sample_weight
            else:
                tree._fit(
                    X,
                    y,
                    sample_weight=sample_weight,
                    check_input=False,
                    missing_values_in_feature_mask=missing_values_in_feature_mask,
                )
                return tree, None

        #@_fit_context(prefer_skip_nested_validation=True)
        def fit(self, X, y, sample_weight=None):
            """
            Build a forest of trees from the training set (X, y).

            Parameters
            ----------
            X : {array-like, sparse matrix} of shape (n_samples, n_features)
                The training input samples. Internally, its dtype will be converted
                to ``dtype=np.float32``. If a sparse matrix is provided, it will be
                converted into a sparse ``csc_matrix``.

            y : array-like of shape (n_samples,) or (n_samples, n_outputs)
                The target values (class labels in classification, real numbers in
                regression).

            sample_weight : array-like of shape (n_samples,), default=None
                Sample weights. If None, then samples are equally weighted. Splits
                that would create child nodes with net zero or negative weight are
                ignored while searching for a split in each node. In the case of
                classification, splits are also ignored if they would result in any
                single class carrying a negative weight in either child node.

            Returns
            -------
            self : object
                Fitted estimator.
            """
            # Validate or convert input data
            if issparse(y):
                raise ValueError("sparse multilabel-indicator for y is not supported.")

            X, y = self._validate_data(
                X,
                y,
                multi_output=True,
                accept_sparse="csc",
                dtype=DTYPE,
                force_all_finite=False,
            )
            # _compute_missing_values_in_feature_mask checks if X has missing values and
            # will raise an error if the underlying tree base estimator can't handle missing
            # values. Only the criterion is required to determine if the tree supports
            # missing values.
            estimator = type(self.estimator)(criterion=self.criterion)
            missing_values_in_feature_mask = (
                estimator._compute_missing_values_in_feature_mask(
                    X, estimator_name=self.__class__.__name__
                )
            )

            if sample_weight is not None:
                sample_weight = _check_sample_weight(sample_weight, X)

            if issparse(X):
                # Pre-sort indices to avoid that each individual tree of the
                # ensemble sorts the indices.
                X.sort_indices()

            y = np.atleast_1d(y)
            if y.ndim == 2 and y.shape[1] == 1:
                warn(
                    (
                        "A column-vector y was passed when a 1d array was"
                        " expected. Please change the shape of y to "
                        "(n_samples,), for example using ravel()."
                    ),
                    DataConversionWarning,
                    stacklevel=2,
                )

            if y.ndim == 1:
                # reshape is necessary to preserve the data contiguity against vs
                # [:, np.newaxis] that does not.
                y = np.reshape(y, (-1, 1))

            if self.criterion == "poisson":
                if np.any(y < 0):
                    raise ValueError(
                        "Some value(s) of y are negative which is "
                        "not allowed for Poisson regression."
                    )
                if np.sum(y) <= 0:
                    raise ValueError(
                        "Sum of y is not strictly positive which "
                        "is necessary for Poisson regression."
                    )

            self._n_samples, self.n_outputs_ = y.shape

            y, expanded_class_weight = self._validate_y_class_weight(y)

            if getattr(y, "dtype", None) != DOUBLE or not y.flags.contiguous:
                y = np.ascontiguousarray(y, dtype=DOUBLE)

            if expanded_class_weight is not None:
                if sample_weight is not None:
                    sample_weight = sample_weight * expanded_class_weight
                else:
                    sample_weight = expanded_class_weight

            if not self.bootstrap and self.max_samples is not None:
                raise ValueError(
                    "`max_sample` cannot be set if `bootstrap=False`. "
                    "Either switch to `bootstrap=True` or set "
                    "`max_sample=None`."
                )
            elif self.bootstrap:
                n_samples_bootstrap = self._get_n_samples_bootstrap(
                    n_samples=X.shape[0], max_samples=self.max_samples
                )
            else:
                n_samples_bootstrap = None

            self._n_samples_bootstrap = n_samples_bootstrap

            self._validate_estimator()

            if not self.bootstrap and self.oob_score:
                raise ValueError("Out of bag estimation only available if bootstrap=True")

            random_state = check_random_state(self.random_state)

            if not self.warm_start or not hasattr(self, "estimators_"):
                # Free allocated memory, if any
                self.estimators_ = []

            n_more_estimators = self.n_estimators - len(self.estimators_)

            if n_more_estimators < 0:
                raise ValueError(
                    "n_estimators=%d must be larger or equal to "
                    "len(estimators_)=%d when warm_start==True"
                    % (self.n_estimators, len(self.estimators_))
                )

            elif n_more_estimators == 0:
                warn(
                    "Warm-start fitting without increasing n_estimators does not "
                    "fit new trees."
                )
            else:
                if self.warm_start and len(self.estimators_) > 0:
                    # We draw from the random state to get the random state we
                    # would have got if we hadn't used a warm_start.
                    random_state.randint(MAX_INT, size=len(self.estimators_))

                trees = [
                    self._make_estimator(append=False, random_state=random_state)
                    for i in range(n_more_estimators)
                ]

                # Parallel loop: we prefer the threading backend as the Cython code
                # for fitting the trees is internally releasing the Python GIL
                # making threading more efficient than multiprocessing in
                # that case. However, for joblib 0.12+ we respect any
                # parallel_backend contexts set at a higher level,
                # since correctness does not rely on using threads.
                trees = Parallel(
                    n_jobs=self.n_jobs,
                    verbose=self.verbose,
                    prefer="threads",
                )(
                    delayed(self._parallel_build_trees)(
                        t,
                        self.bootstrap,
                        X,
                        y,
                        sample_weight,
                        i,
                        len(trees),
                        verbose=self.verbose,
                        class_weight=self.class_weight,
                        n_samples_bootstrap=n_samples_bootstrap,
                        missing_values_in_feature_mask=missing_values_in_feature_mask,
                    )
                    for i, t in enumerate(trees)
                )

                tree_models, index_counts = map(list, zip(*trees))
                # Collect newly grown trees
                self.estimators_.extend(tree_models)
                self.index_counts = index_counts

            if self.oob_score and (
                n_more_estimators > 0 or not hasattr(self, "oob_score_")
            ):
                y_type = type_of_target(y)
                if y_type == "unknown" or (
                    self._estimator_type == "classifier"
                    and y_type == "multiclass-multioutput"
                ):
                    # FIXME: we could consider to support multiclass-multioutput if
                    # we introduce or reuse a constructor parameter (e.g.
                    # oob_score) allowing our user to pass a callable defining the
                    # scoring strategy on OOB sample.
                    raise ValueError(
                        "The type of target cannot be used to compute OOB "
                        f"estimates. Got {y_type} while only the following are "
                        "supported: continuous, continuous-multioutput, binary, "
                        "multiclass, multilabel-indicator."
                    )

                if callable(self.oob_score):
                    self._set_oob_score_and_attributes(
                        X, y, scoring_function=self.oob_score
                    )
                else:
                    self._set_oob_score_and_attributes(X, y)

            # Decapsulate classes_ attributes
            if hasattr(self, "classes_") and self.n_outputs_ == 1:
                self.n_classes_ = self.n_classes_[0]
                self.classes_ = self.classes_[0]

            return self
