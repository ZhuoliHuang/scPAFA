import pandas as pd
import anndata as ad
import numpy as np
from typing import Dict, List


def map_genes_to_positions(gene_dict, gene_array):
    # a function to map gene in dict to a index of gene array
    gene_dict = gene_dict.copy()
    # 创建一个反向映射，将基因名字映射到在数组中的位置
    gene_positions = {gene: index for index, gene in enumerate(gene_array)}

    # 更新字典的值为基因在数组中的位置
    for key, value in gene_dict.items():
        if isinstance(value, list):
            # 如果值是一个基因列表，映射所有基因到位置
            gene_dict[key] = [gene_positions[g] for g in value]
        else:
            # 如果值是单个基因，映射它到位置
            gene_dict[key] = gene_positions.get(value, None)

    return gene_dict


def filter_dict_by_intersection(input_dict: Dict[str, List[str]], 
                                target_list: List[str],
                                min_gene_num: int) -> Dict[str, List[str]]:
    """
    Filter a dictionary of pathways based on the minimum overlap with a target gene list.

    Parameters:
        input_dict (dict): A dictionary where keys are pathway names and values are lists of genes.
        target_list (list): A list of target genes to compare against pathway genes.
        min_gene_num (int): The minimum number of overlapping genes required to keep a pathway.

    Returns:
        dict: A filtered dictionary containing pathways with sufficient overlap with the target gene list.
    """
    filtered_dict = {}
    filtered_out_count = 0  # a_counter_for_filtered_cell

    for key, value in input_dict.items():
        intersection_count = len(set(value) & set(target_list))
        
        if intersection_count >= min_gene_num:
            filtered_dict[key] = value
        else:
            filtered_out_count += 1
    print(f"Filtered out {filtered_out_count} pathways")
    return filtered_dict

#function_for_pathway_input_for_UCell
def generate_pathway_input(adata: ad.AnnData,
                          pathway_dict: Dict[str, List[str]],
                          min_overlap_gene: int = 3) -> dict:
    """
    Generate pathway input data for fast_UCell or fast_score_genes analysis.

    Parameters:
        adata (ad.AnnData): An AnnData object containing gene expression data.
        pathway_dict (Dict[str, List[str]]): A dictionary where keys are pathway names and values are lists of genes.
        min_overlap_gene (int, optional): The minimum overlap of genes required to include a pathway. Default is 3.

    Returns:
        dict: A dictionary containing 3 dicts and 2 1-D arrays.
        'pathway_dict_length' : a dict of original length of pathways(fast_UCell)
        'pathway_dict_filtered_gene': a dict of genes overlap in adata.var_names and pathways(fast_score_genes)
        'intersectgene': a array of union set of all gene in 'pathway_dict_filtered_gene'
        'intersect_position':the index of 'intersectgene' in adata.var_names (fast_UCell)
        pathway_dict_filtered_position':a dict map gene in 'pathway_dict_filtered_gene' to index of 'intersectgene' (fast_UCell)
    """
    if not isinstance(adata, ad.AnnData):
        raise TypeError("Input 'adata' must be an AnnData object.")
    if not isinstance(pathway_dict, dict):
        raise TypeError("Input 'pathway_dict' must be a dictionary.")
    if not isinstance(min_overlap_gene, int):
        raise TypeError("Input 'min_overlap_gene' must be an integer.")
    
    #直接过滤掉通路，不影响通路的基因
    pathway_dict_filtered = filter_dict_by_intersection(pathway_dict,target_list=list(adata.var_names),min_gene_num = min_overlap_gene)
    #记录原始通路的基因数
    pathway_dict_length = {key: len(value) for key, value in pathway_dict_filtered.items()}
    
    #将通路dict展开
    pathway_genes = [item for value in pathway_dict_filtered.values() for item in (value if isinstance(value, list) else [value])]
    #去除重复的基因
    pathway_genes = np.array(list(set(pathway_genes)),dtype = 'object')
    #筛选与adata有overlap的通路基因
    pathway_genes = pathway_genes[pd.Series(pathway_genes).isin(adata.var_names)]
    
    #创建一个每个通路都只剩下和adata的基因交集的字典(用于addmoudle_score)
    filtered_gene_dict = {key: [gene for gene in value if gene in pathway_genes] for key, value in pathway_dict_filtered.items()}
    
    #找到交集基因在adata.X中的列位置，并且将交集基因按照从左到右排序
    pathway_genes_position = np.where(adata.var_names.isin(pathway_genes))[0]
    pathway_genes = np.array(adata.var_names[pathway_genes_position])
    
    #将交集基因字典中的基因转成在pathway_genes中的排序index
    filtered_gene_dict_position = map_genes_to_positions(filtered_gene_dict,pathway_genes)
    
    # a_dict_for_result
    result_dict = {
        'pathway_dict_length': pathway_dict_length,
        'pathway_dict_filtered_gene': filtered_gene_dict,
        'pathway_dict_filtered_position':filtered_gene_dict_position,
        'intersectgene': pathway_genes,
        'intersect_position':pathway_genes_position
    }
    
    print(str(len(pathway_dict_filtered.keys()))+' pathways passed QC')
    print('The maxRank must >= '+str(max(result_dict['pathway_dict_length'].values())),'(The genes number of the longest pathway)')
    return result_dict