import os
import numpy as np
import pandas as pd
#for_adata
import anndata as ad
import scipy
#multi_processing_module
import pathos
from tqdm import tqdm
from math import ceil

#
# UCell:step 1 rank, step 2 score 
#

#UCell:step 1 rank
def fast_UCell_rank(
            adata:ad.AnnData = None,
            maxRank:int = 1500,
            n_cores_rank:int = None,
            rank_batch_size:int = 100000):
    
    """
    Perform a fast UCell analysis on single-cell RNA-seq data. Step1 ranking.
    For each cell, rank gene expression from high to low.

    Parameters
    ----------
    adata : ad.AnnData,Required
        An AnnData object containing the single-cell RNA-seq data.
    maxRank : int, optional
        The maximum rank for UCell. Default is 1500. maxRank must bigger than the gene numbers of the longest pathway.
    n_cores_rank : int, optional
        The number of CPU cores to use for ranking.
        Please note that adjusting this parameter should take into account rank_batch_size.
        rank_batch_size/n_cores_rank above 10,000 is recommanded. 
        Having more CPU cores does not necessarily guarantee performance improvement.
    rank_batch_size : int, optional
        The batch size for ranking computation. Default is 100000.
        A smaller batch will take less memory.
        Will extract a batch of rank_batch_size cells, split them into smaller chunks for multiprocessing.
        
    Returns
    -------
    A rank scipy sparse csr matrix with row as cells and columns as genes.

    Raises
    ------
    ValueError
        If the provided adata is not a valid AnnData object.
        If maxRank is not a positive integer,or it is smaller than the max gene number among pathways.
        If n_cores_rank is not a positive integer.
        If rank_batch_size is not a positive integer.
    """
    
    if not isinstance(adata, ad.AnnData):
        raise ValueError("adata must be a valid AnnData object.")
    
    if not isinstance(n_cores_rank, int) or n_cores_rank <= 0 or n_cores_rank > pathos.multiprocessing.cpu_count():
        raise ValueError("n_cores_rank must be a positive integer and smaller than max cpu_count")

    if not isinstance(rank_batch_size, int) or rank_batch_size <= 0:
        raise ValueError("rank_batch_size must be a positive integer.")
        
    def worker(i):
        #chunk 
        start = i * rank_chunk_size
        end = min((i + 1) * rank_chunk_size, len(count_batch))
        count_df = count_batch[start:end]
        #rank
        rank_array = scipy.stats.rankdata((count_df*(-1)), axis=1)
        rank_array[rank_array > maxRank] = 0
        rank_sparse_matrix = scipy.sparse.csr_matrix(rank_array)
        return rank_sparse_matrix

    count = pd.DataFrame.sparse.from_spmatrix(adata.X)
    num_batches = len(count) //  rank_batch_size
    if len(count) %  rank_batch_size != 0:
        num_batches += 1
    
    num_chunks = n_cores_rank
    rank_chunk_size = ceil(rank_batch_size/num_chunks)
    
    print('Step1 generate rank matrix')
    print(str(num_batches) + ' batches need to rank, with each max '+str(rank_batch_size)+' cells')
    final_csr_matrix_list=[]

    for i_batch in range(num_batches):
        print('Processing_batch_'+str((i_batch+1)))
        start = i_batch * rank_batch_size
        end = min((i_batch + 1) * rank_batch_size,len(count))
        count_batch = count[start:end]
    
        csr_matrix_list = []
    
        with pathos.pools.ProcessPool(nodes=n_cores_rank) as pool:
            results = list(tqdm(pool.imap(worker, range(num_chunks), chunksize=1), total=num_chunks, desc="Ranking Chunks"))
    
        #results.sort(key=lambda x: x[0])
        # Extract rank_sparse_matrix from sorted results
        
        csr_matrix_list = [result for result in results]
        csr_matrix_batch = scipy.sparse.vstack(csr_matrix_list)
        final_csr_matrix_list.append(csr_matrix_batch)
    
    final_csr_matrix = scipy.sparse.vstack(final_csr_matrix_list)
    
    print('Rank done')
    print('The output rank matrix can be used to calculate UCell score on different pathways set')
    print('The maxRank parameter use with this rank matrix in fast_UCell_score() must <= '+str(maxRank))
    return(final_csr_matrix)


#UCell:step 2 score
def fast_UCell_score(
            cell_index = None,
            rankmatrix:scipy.sparse.csr_matrix = None,
            maxRank:int = None,
            input_dict:dict = None,
            score_batch_size:int = 100000,
            n_cores_score:int = None):
    
    """
    Perform a fast UCell analysis on single-cell RNA-seq data.

    Parameters
    ----------
    cell_index : the cell index of adata used in fast_UCell_rank() for example: adata.obs.index
        An array like
    maxRank : int, optional
        The maximum rank for UCell,must <= maxRank used in fast_UCell_rank().
    input_dict : dict, Required
        A dictionary from generate_pathway_input function, containing input data.
    score_batch_size : int, optional
        The batch size for scoring. Default is 100000.
        Will extract a batch of score_batch_size cells.The pathways are distributed across various cores, and calculations are performed on these cells.
    n_cores_score : int
        The number of CPU cores to use for scoring.
        For hundreds of pathways, 4~8 cores are fast enough.
        For scenarios with a high number of pathways, increasing the number of CPU cores will significantly reduce the runtime.
    Returns
    -------
    A UCell score pandas dataframe with row index as cells and columns as pathways.

    Raises
    ------
    ValueError
        If input_dict is not a dict.
        If maxRank is not a positive integer,or it is smaller than the max gene number among pathways.
        If n_cores_score is not a positive integer.
    """   
    if not isinstance(maxRank, int) or maxRank <= 0:
        raise ValueError("maxRank must be a positive integer.")
        
    if not isinstance(input_dict, Dict):
        raise ValueError("input_dict must be a dict.")
    
    if maxRank < max(input_dict['pathway_dict_length'].values()):
        raise ValueError("maxRank is smaller than "+str(max(test_input['pathway_dict_length'].values()))+
                                                        ' which is the genes number of the longest pathway')
        
    if not isinstance(score_batch_size, int) or score_batch_size <= 0:
        raise ValueError("score_batch_size must be a positive integer.")
    
    if not isinstance(n_cores_score, int) or n_cores_score <= 0 or n_cores_score > pathos.multiprocessing.cpu_count():
        raise ValueError("n_cores_score must be a positive integer and smaller than max cpu_count")

    
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
    keys = list(pathway_pos.keys())
    
    num_batches = final_csr_matrix.shape[0] // score_batch_size
    if final_csr_matrix.shape[0] %  score_batch_size != 0:
        num_batches += 1

    print('step2 calculating Score')
    print(str(num_batches) + ' batches need to score, with each max '+str(score_batch_size)+' cells')
    
    final_score_list=[]
    
    for i_batch in range(num_batches):
        print('processing_batch_'+str((i_batch+1)))
        start = i_batch * score_batch_size
        end = min((i_batch + 1) * score_batch_size,final_csr_matrix.shape[0])
        final_csr_matrix_batch = final_csr_matrix[start:end,:]
        
        with pathos.pools.ProcessPool(nodes=n_cores_score) as pool:
            score_results = list(tqdm(pool.imap(process_key,keys,chunksize=ceil(len(keys)/n_cores_score)), total=len(keys), desc="Pathways"))
        Ucell_dataframe = pd.concat(score_results,axis=1)
        final_score_list.append(Ucell_dataframe)
        
    print('Ucell_done')
    print('Outputing_dataframe')
    UCell_score_df = pd.concat(final_score_list, axis=0, ignore_index=True)
    UCell_score_df.index = cell_index
    #Ucell_dataframe = pd.concat(score_results,axis=1)
    return UCell_score_df