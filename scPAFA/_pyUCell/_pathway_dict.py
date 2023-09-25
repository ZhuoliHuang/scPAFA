import pandas as pd
import anndata as ad
import numpy as np
from typing import Dict, List

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
    Generate pathway input data for UCell analysis.

    Parameters:
        adata (ad.AnnData): An AnnData object containing gene expression data.
        pathway_dict (Dict[str, List[str]]): A dictionary where keys are pathway names and values are lists of genes.
        min_overlap_gene (int, optional): The minimum overlap of genes required to include a pathway. Default is 3.

    Returns:
        dict: A dictionary containing pathway dataframes, pathway length information, and intersecting genes.
    """
    if not isinstance(adata, ad.AnnData):
        raise TypeError("Input 'adata' must be an AnnData object.")
    if not isinstance(pathway_dict, dict):
        raise TypeError("Input 'pathway_dict' must be a dictionary.")
    if not isinstance(min_overlap_gene, int):
        raise TypeError("Input 'min_overlap_gene' must be an integer.")
    
    pathway_dict_filtered = filter_dict_by_intersection(pathway_dict,target_list=list(adata.var_names),min_gene_num = min_overlap_gene)
    Dict_dataframe = pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in pathway_dict_filtered.items() ]))
    info_Dict_dataframe = pd.DataFrame({'length':Dict_dataframe.count()})
    Dict_dataframe[-(Dict_dataframe.isin(adata.var_names))] = np.nan
    info_Dict_dataframe['length_inter'] = Dict_dataframe.count()
    gene_name = Dict_dataframe.stack()
    #intersect_gene_of_all_pathway_and_adata.var
    intersectgene = gene_name.drop_duplicates().values

    # a_dict_for_result
    result_dict = {
        'pathway_dataframe': Dict_dataframe,
        'pathway_length_dataframe': info_Dict_dataframe,
        'intersectgene': intersectgene
    }
    
    return result_dict