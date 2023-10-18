import numpy as np
import pandas as pd
import anndata as ad
import pathos
import scipy
from scipy.sparse import issparse
from typing import Sequence, Optional,Union
from tqdm import tqdm
from numpy import random
from math import ceil

def _sparse_nanmean(X, axis):
    """
    np.nanmean equivalent for sparse matrices
    """
    if not issparse(X):
        raise TypeError("X must be a sparse matrix")

    # count the number of nan elements per row/column (dep. on axis)
    Z = X.copy()
    Z.data = np.isnan(Z.data)
    Z.eliminate_zeros()
    n_elements = Z.shape[axis] - Z.sum(axis)

    # set the nans to 0, so that a normal .sum() works
    Y = X.copy()
    Y.data[np.isnan(Y.data)] = 0
    Y.eliminate_zeros()

    # the average
    s = Y.sum(axis, dtype='float64')  # float64 for score_genes function compatibility)
    m = s / n_elements

    return m

def _check_use_raw(adata: ad.AnnData, use_raw: Union[None, bool]) -> bool:
    """
    Normalize checking `use_raw`.

    My intentention here is to also provide a single place to throw a deprecation warning from in future.
    """
    if use_raw is not None:
        return use_raw
    else:
        if adata.raw is not None:
            return True
        else:
            return False

def fast_sctl_score(
    adata:ad.AnnData = None,
    input_dict:dict = None,
    n_cores_score:int = None,
    score_batch_size:int = 100000,
    random_seed: int = 0,
    ctrl_size: int = 50,
    n_bins: int = 25,
    use_raw: Optional[bool] = None,
    gene_pool: Optional[Sequence[str]] = None):
    
    """
    Fast multi-core version of scanpy.tl.score_genes function.
    The score is the average expression of a set of genes subtracted with the average expression of a reference set of genes. 
    The reference gene set is randomly sampled from the gene_pool for each binned expression value.
    
    https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.score_genes.html#scanpy.tl.score_genes
    Satija et al. (2015), Spatial reconstruction of single-cell gene expression data.
    
    Parameters
    ----------
    adata : ad.AnnData,Required
        An AnnData object containing the single-cell RNA-seq data.
    input_dict : dict, Required
        A dictionary from generate_pathway_input function, containing input data.
    n_cores_score : int,Required
        The number of CPU cores to use for scoring.
        For scenarios with a high number of pathways, increasing the number of CPU cores will significantly reduce the runtime.
    score_batch_size:int,Required
        The batch size for scoring.Default is 100000.
        After calculate binned expression of all genes(This step can not be sliced to run)
        Will extract a batch of score_batch_size cells.
        The pathways are distributed across various cores, and calculations are performed on these cells.
    random_seed:int
        seed to replicate result.Default is 0.
        This is used in randomly selecting reference gene.
    ctrl_size:int
        Number of reference genes to be sampled from each bin. Default is 50.
    n_bins:int 
        Number of expression level bins for sampling.Default is 25.
    use_raw:Optional[bool]
        Whether to use raw attribute of adata. Defaults to True if .raw is present.
    gene_pool: Optional[Sequence[str]]
        Genes for sampling the reference set. Default is all genes.
    
    Returns
    -------
    A rank scipy sparse csr matrix with row as cells and columns as genes.
    """
    # Check if adata is a valid AnnData object
    if not isinstance(adata, ad.AnnData):
        raise ValueError("adata must be a valid AnnData object.")

    if not isinstance(input_dict, Dict):
        raise ValueError("input_dict must be a dict.")
    
    if not isinstance(score_batch_size, int) or score_batch_size <= 0:
        raise ValueError("score_batch_size must be a positive integer.")
    
    # Check if n_cores_score is positive integer
    if not isinstance(n_cores_score, int) or n_cores_score <= 0 or n_cores_score > pathos.multiprocessing.cpu_count():
        raise ValueError("n_cores_score must be a positive integer and smaller than max cpu_count")
    
    
    def select_control_genes(pathway_genelist):
        np.random.seed(random_seed)
        gene_list = pathway_genelist
        gene_list = set(gene_list)
        control_genes = set()
        # now pick `ctrl_size` genes from every cut
        for cut in np.unique(obs_cut.loc[list(gene_list)]):
            r_genes = np.array(obs_cut[obs_cut == cut].index)
            np.random.shuffle(r_genes)
            # uses full r_genes if ctrl_size > len(r_genes)
            control_genes.update(set(r_genes[:ctrl_size]))

        control_genes = control_genes - gene_list
        control_genes = list(control_genes)
        return control_genes    
    
    def calculate_score_for_pathway(pathway_key):
        query_genes = np.where(all_gene.isin(pathway_dict_filtered_gene[pathway_key]))[0]
        control_genes = select_control_genes(pathway_genelist=pathway_dict_filtered_gene[pathway_key])
        control_genes = np.where(all_gene.isin(control_genes))[0]

        X_list = X_to_use_batch[:, query_genes]
        X_control = X_to_use_batch[:, control_genes]

        if sparse_or_not:
            X_list = np.array(_sparse_nanmean(X_list, axis=1)).flatten()
            X_control = np.array(_sparse_nanmean(X_control, axis=1)).flatten()
        else:
            X_list = np.nanmean(X_list, axis=1, dtype='float64')
            X_control = np.nanmean(X_control, axis=1, dtype='float64')

        score = X_list - X_control
        score = pd.Series(score, dtype='float64', name=pathway_key)

        return score
    
    
    use_raw = _check_use_raw(adata, use_raw)
    var_names = adata.raw.var_names if use_raw else adata.var_names

    if gene_pool is None:
        gene_pool = list(var_names)
    else:
        gene_pool = [x for x in gene_pool if x in var_names]
    if not gene_pool:
        raise ValueError("No valid genes were passed for reference set.")
    
    _adata = adata.raw if use_raw else adata
    _adata_subset = (_adata[:, gene_pool] if len(gene_pool) < len(_adata.var_names) else _adata)
    
    sparse_or_not = issparse(_adata_subset.X)
    
    if sparse_or_not:
        obs_avg = pd.Series(
            np.array(_sparse_nanmean(_adata_subset.X, axis=0)).flatten(),
            index=gene_pool,
        )  # average expression of genes
    else:
        obs_avg = pd.Series(
            np.nanmean(_adata_subset.X, axis=0), index=gene_pool
        )  # average expression of genes

    obs_avg = obs_avg[np.isfinite(obs_avg)]
    n_items = int(np.round(len(obs_avg) / (n_bins - 1)))
    obs_cut = obs_avg.rank(method='min') // n_items
    
    X_to_use = _adata_subset.X
    all_gene = _adata_subset.var_names
    pathway_dict_filtered_gene = input_dict['pathway_dict_filtered_gene']
    
    num_batches = X_to_use.shape[0] // score_batch_size
    if X_to_use.shape[0] %  score_batch_size != 0:
        num_batches += 1
    
    print('calculating Score')
    print(str(num_batches) + ' batches need to score, with each max '+str(score_batch_size)+' cells')
    
    keys = list(pathway_dict_filtered_gene.keys())
    final_score_list=[]
    
    for i_batch in range(num_batches):
        print('processing_batch_'+str((i_batch+1)))
        start = i_batch * score_batch_size
        end = min((i_batch + 1) * score_batch_size,X_to_use.shape[0])
        X_to_use_batch = X_to_use[start:end,:]
        with pathos.pools.ProcessPool(nodes=n_cores_score) as pool:
            score_results = pool.map(calculate_score_for_pathway,keys)
        score_dataframe = pd.concat(score_results,axis=1)
        final_score_list.append(score_dataframe)
        
    print('score_done')
    print('Outputing_dataframe')
    score_df = pd.concat(final_score_list, axis=0, ignore_index=True)
    score_df.index = _adata_subset.obs.index
    
    return score_df