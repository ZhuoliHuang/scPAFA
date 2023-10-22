import numpy as np
import pandas as pd
import anndata as ad
import scipy
import pathos
from tqdm import tqdm
from math import ceil
from scipy.sparse import issparse
from typing import Sequence

# UCell:step 1 rank, step 2 score 

#UCell:step 1 rank
def fast_ucell_rank(
            adata:ad.AnnData,
            n_cores_rank:int,
            maxRank:int = 1500,
            rank_batch_size:int = 100000):
    
    """
    Perform a fast UCell analysis on single-cell RNA-seq data. Step1 ranking.
    For each cell, rank gene expression(from count or lognormoalized data) from high to low.
    UCell algorithm: Andreatta, M., & Carmona, S. J. (2021). UCell: Robust and scalable single-cell gene signature scoring. 

    Parameters
    ----------
    adata : ad.AnnData
        An AnnData object containing the single-cell RNA-seq data, with adata.X as raw count(or log-normalized data) is recommanded. 
        The results are identical for both raw count and log-normalized data. Because scanpy.pp.normalize_total() and scanpy.pp.log1p() do not alter the expression ranking.
    n_cores_rank : int
        The number of CPU cores to use for ranking.
        Please note that adjusting this parameter should take into account rank_batch_size.
        When using raw count, rank_batch_size/n_cores_rank above 10,000 is recommanded. 4~8 cores when rank_batch_size is 100,000 to 200,000 works good.
        Having more CPU cores does not necessarily guarantee performance improvement.
    maxRank : int, optional
        The maximum rank for UCell. Default is 1500. maxRank must bigger than the gene numbers of the longest pathway.
    rank_batch_size : int, optional
        The batch size for ranking computation. Default is 100000.
        A smaller batch will take less memory.
        Using scaled data will cost more memory.
        Will extract a batch of rank_batch_size cells, split them into smaller chunks for multiprocessing.
        
    Returns
    -------
    A rank scipy sparse csr matrix with row as cells and columns as genes.

    Raises
    ------
    ValueError
        If the provided adata is not a valid AnnData object.
        If n_cores_rank is not a positive integer.
        If maxRank is not a positive integer.
        If rank_batch_size is not a positive integer.
    """
    
    if not isinstance(adata, ad.AnnData):
        raise ValueError("adata must be a valid AnnData object.")
    
    if not isinstance(n_cores_rank, int) or n_cores_rank <= 0 or n_cores_rank > pathos.multiprocessing.cpu_count():
        raise ValueError("n_cores_rank must be a positive integer and <= max cpu_count")
    
    if not isinstance(maxRank, int) or maxRank <= 0:
        raise ValueError("maxRank must be a positive integer.")    
    
    if not isinstance(rank_batch_size, int) or rank_batch_size <= 0:
        raise ValueError("rank_batch_size must be a positive integer.")
    
    # a rank function based on scipy.stats.rankdata that will run on each core
    def worker(i):
        #chunk 
        start = i * rank_chunk_size
        end = min((i + 1) * rank_chunk_size, len(count_batch))
        count_df = count_batch[start:end]

        #rank
        rank_array = scipy.stats.rankdata((count_df*(-1)), axis=1)
        #gene expression rank bigger than maxRank will be set to 0 to save memory
        rank_array[rank_array > maxRank] = 0
        rank_sparse_matrix = scipy.sparse.csr_matrix(rank_array)
        return rank_sparse_matrix

    # check the type of adata.X
    if issparse(adata.X):
        count = pd.DataFrame.sparse.from_spmatrix(adata.X)
        print('adata.X is csr sparse matrix')
    else:
        count = adata.X.copy()
        print('adata.X is numpy.ndarray')

    num_batches = len(count) //  rank_batch_size
    if len(count) %  rank_batch_size != 0:
        num_batches += 1
    
    #manually divide by the number of cores
    num_chunks = n_cores_rank
    rank_chunk_size = ceil(rank_batch_size/num_chunks)
    
    print('Step1 generate rank matrix')
    print(str(num_batches) + ' batches need to rank, with each max '+str(rank_batch_size)+' cells')
    final_csr_matrix_list=[]

    for i_batch in range(num_batches):
        
        # divide into batch
        print('Processing_batch_'+str((i_batch+1)))
        start = i_batch * rank_batch_size
        end = min((i_batch + 1) * rank_batch_size,len(count))
        count_batch = count[start:end]
    
        csr_matrix_list = []
        
        #multiprocessing
        with pathos.pools.ProcessPool(nodes=n_cores_rank) as pool:
            results = list(tqdm(pool.imap(worker, range(num_chunks), chunksize=1), total=num_chunks, desc="Ranking Chunks"))

        #merge result from each chunk
        csr_matrix_list = [result for result in results]
        csr_matrix_batch = scipy.sparse.vstack(csr_matrix_list)
        final_csr_matrix_list.append(csr_matrix_batch)
    
    #merge result from each batch
    final_csr_matrix = scipy.sparse.vstack(final_csr_matrix_list)
    
    print('Rank done')
    print('The output rank matrix can be used to calculate UCell score on different pathways set.')
    print('The maxRank parameter use with this rank matrix in fast_ucell_score() must <= '+str(maxRank))
    return(final_csr_matrix)


#UCell:step 2 score
def fast_ucell_score(
            cell_index:Sequence,
            rankmatrix:scipy.sparse.csr_matrix,
            n_cores_score:int,
            input_dict:dict ,
            maxRank:int = 1500,
            score_batch_size:int = 100000):
    
    """
    Perform a fast UCell analysis on single-cell RNA-seq data.
    UCell algorithm: Andreatta, M., & Carmona, S. J. (2021). UCell: Robust and scalable single-cell gene signature scoring.

    Parameters
    ----------
    cell_index : Sequence[str]
        The cell index of adata used in fast_ucell_rank(), for example: list(adata.obs.index).
    rankmatrix : scipy.sparse.csr_matrix
        The result of fast_ucell_rank().
    n_cores_score : int
        The number of CPU cores to use for scoring.
        For hundreds of pathways, 4~8 cores are fast enough.
        For scenarios with a high number of pathways, increasing the number of CPU cores will significantly reduce the runtime.
    input_dict : dict
        A dictionary from generate_pathway_input(), containing input data.
    maxRank : int, optional
        The maximum rank for UCell,must <= maxRank used in fast_ucell_rank().Default is 1500.
    score_batch_size : int, optional
        The batch size for scoring. Default is 100000.
        Will extract a batch of score_batch_size cells.The pathways are distributed across various cores, and calculations are performed on these cells.
    Returns
    -------
    A UCell score pandas dataframe with row index as cells and columns as pathways.

    Raises
    ------
    ValueError:
        If `cell_index` is not an array-like of strings.
        If `rankmatrix` is not a `scipy.sparse.csr_matrix`.
        If `n_cores_score` is not a positive integer or exceeds the maximum CPU count.
        If `input_dict` is not a dictionary.
        If `maxRank` is smaller than the maximum gene count among pathways or not a positive integer.
        If `score_batch_size` is not a positive integer.
    """   
    if not isinstance(cell_index, Sequence):
        raise ValueError("cell_index must be a list of strings.")
    
    if not isinstance(rankmatrix, scipy.sparse.csr_matrix):
        raise ValueError("rankmatrix must be scipy.sparse.csr_matrix.")

    if not len(cell_index) == rankmatrix.shape[0]:
        raise ValueError("The length of cell_index must match the rows number of rankmatrix.")

    if not isinstance(n_cores_score, int) or n_cores_score <= 0 or n_cores_score > pathos.multiprocessing.cpu_count():
        raise ValueError("n_cores_score must be a positive integer and smaller than max cpu_count.")

    if not isinstance(input_dict, dict):
        raise ValueError("input_dict must be a dict.")
    
    if maxRank < max(input_dict['pathway_dict_length'].values()):
        raise ValueError("maxRank is smaller than "+str(max(input_dict['pathway_dict_length'].values()))+
                                                        ' which is the genes number of the longest pathway.')
    if not isinstance(maxRank, int) or maxRank <= 0:
        raise ValueError("maxRank must be a positive integer.")

    if not isinstance(score_batch_size, int) or score_batch_size <= 0:
        raise ValueError("score_batch_size must be a positive integer.")
    
    #mission run for each pathway
    def process_key(key):
        #print('processing_'+ key)
        length = pathway_length[key]
    
        num1 = (length*(length+1))/2
        num2 = maxRank*length
    
        position = np.array(pathway_pos[key])
        subset_final_csr_matrix = final_csr_matrix_batch[:,position]
    
        special_index = np.array(subset_final_csr_matrix.sum(axis = 1)).ravel() != 0
        overmaxRank = length - subset_final_csr_matrix.getnnz(axis = 1)

        final_result = np.array(subset_final_csr_matrix.sum(axis = 1)).ravel()
        final_result[special_index] = final_result[special_index]+(maxRank+1)*overmaxRank[special_index]
        final_result[special_index] = 1 - (final_result[special_index] - num1) / num2
        final_result = pd.Series(final_result,name=key)
        return final_result
    
    print('Subset rank matrix by overlap genes in pathways')
    final_csr_matrix = rankmatrix.copy()
    final_csr_matrix = final_csr_matrix[:, input_dict['intersect_position']]

    print('Rank above maxRank to 0')
    final_csr_matrix.data[final_csr_matrix.data > maxRank] = 0
    final_csr_matrix.eliminate_zeros()
    
    pathway_pos = input_dict['pathway_dict_filtered_position']
    pathway_length = input_dict['pathway_dict_length']
    
    #get all pathway names
    keys = list(pathway_pos.keys())
    
    # divide into batch
    num_batches = final_csr_matrix.shape[0] // score_batch_size
    if final_csr_matrix.shape[0] %  score_batch_size != 0:
        num_batches += 1

    print('step2 calculating Score')
    print(str(num_batches) + ' batches need to score, with each max '+str(score_batch_size)+' cells')
    
    final_score_list=[]
    
    for i_batch in range(num_batches):
        # divide into batch
        print('processing_batch_'+str((i_batch+1)))
        start = i_batch * score_batch_size
        end = min((i_batch + 1) * score_batch_size,final_csr_matrix.shape[0])
        final_csr_matrix_batch = final_csr_matrix[start:end,:]
        
        #multi_processing
        with pathos.pools.ProcessPool(nodes=n_cores_score) as pool:
            score_results = list(tqdm(pool.imap(process_key,keys,chunksize=ceil(len(keys)/n_cores_score)), total=len(keys), desc="Pathways"))

        #merge result for each batch
        Ucell_dataframe = pd.concat(score_results,axis=1)
        final_score_list.append(Ucell_dataframe)
        
    print('Ucell_done')
    print('Outputing_dataframe')
    #merge all result
    UCell_score_df = pd.concat(final_score_list, axis=0, ignore_index=True)
    UCell_score_df.index = cell_index

    return UCell_score_df