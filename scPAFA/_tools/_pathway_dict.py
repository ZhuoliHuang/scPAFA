import pandas as pd
import anndata as ad
import numpy as np
from typing import Dict, List


def map_genes_to_positions(gene_dict, gene_array):
    # a function to map gene in dict to a index of gene array
    gene_dict = gene_dict.copy()
    # Create a reverse mapping that associates gene names with their positions in the array
    gene_positions = {gene: index for index, gene in enumerate(gene_array)}

    # Update the dictionary's values to be the positions of genes in the array
    for key, value in gene_dict.items():
        if isinstance(value, list):
             # If the value is a list of genes, map all genes to their positions
            gene_dict[key] = [gene_positions[g] for g in value]
        else:
            # If the value is a single gene, map it to its position
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

#function_for_pathway_input_for_scoring
def generate_pathway_input(adata: ad.AnnData,
                          pathway_dict: Dict[str, List[str]],
                          min_overlap_gene: int = 3) -> dict:
    """
    Generate pathway input data for fast_ucell or fast_sctl_score analysis.
    Parameters
    ----------
    adata : ad.AnnData
        An AnnData object containing gene expression data.
    pathway_dict : Dict[str, List[str]]
        A dictionary where keys are pathway names and values are lists of genes.
    min_overlap_gene : int, optional
        The minimum overlap of genes required to include a pathway. Default is 3.
    Returns
    -------
    dict
        A dictionary containing the following keys and values:
        - 'pathway_dict_length' (dict): Original length of pathways (for fast_ucell).
        - 'pathway_dict_filtered_gene' (dict): Genes that overlap in adata.var_names and pathways (for fast_score_genes).
        - 'intersectgene' (array): Union set of all genes in 'pathway_dict_filtered_gene'.
        - 'intersect_position' (array): Index of 'intersectgene' in adata.var_names (for fast_UCell).
        - 'pathway_dict_filtered_position' (dict): Mapping of genes in 'pathway_dict_filtered_gene' to the index of 'intersectgene' (for fast_UCell).
    """
    if not isinstance(adata, ad.AnnData):
        raise TypeError("Input 'adata' must be an AnnData object.")
    if not isinstance(pathway_dict, dict):
        raise TypeError("Input 'pathway_dict' must be a dictionary.")
    if not isinstance(min_overlap_gene, int):
        raise TypeError("Input 'min_overlap_gene' must be an integer.")
    
    # Step1: Filter out pathways directly without affecting genes in the pathway
    pathway_dict_filtered = filter_dict_by_intersection(pathway_dict,target_list=list(adata.var_names),min_gene_num = min_overlap_gene)
    
    # Step2: Record the number of genes in the original pathways (useful for fast_ucell)
    pathway_dict_length = {key: len(value) for key, value in pathway_dict_filtered.items()}
    
    # Step3: Flatten the pathway dictionary and Remove duplicate genes
    pathway_genes = [item for value in pathway_dict_filtered.values() for item in (value if isinstance(value, list) else [value])]
    pathway_genes = np.array(list(set(pathway_genes)),dtype = 'object')

    # Step4: Select pathway genes that overlap with 'adata.var_names'
    pathway_genes = pathway_genes[pd.Series(pathway_genes).isin(adata.var_names)]
    
    # Step5: Create a dictionary for each pathway with only genes that intersect with 'adata.var_names' (for fast_sctl_score)
    filtered_gene_dict = {key: [gene for gene in value if gene in pathway_genes] for key, value in pathway_dict_filtered.items()}
    
    # Step6: Find the positions of intersecting genes in 'adata.X' and sort them from left to right
    pathway_genes_position = np.where(adata.var_names.isin(pathway_genes))[0]
    pathway_genes = np.array(adata.var_names[pathway_genes_position])
    
    # Step7: Convert genes in the intersecting gene dictionary to their sorted index in 'pathway_genes'
    filtered_gene_dict_position = map_genes_to_positions(filtered_gene_dict,pathway_genes)
    
    # Step7: Output a dict for result
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