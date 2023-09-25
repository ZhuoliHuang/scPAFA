import numpy as np
import pandas as pd

# A QC funtion for selecting usable sample and view(celltype),based on cell number and sample number

def pseudobulk_qc(
    metadata:pd.DataFrame = None,
    min_cell_number_per_sample: int = 10,
    min_sample_per_view: int = 15,
    sample_column:str = None,
    view_column:str = None):

    """
    Perform quality control for selecting usable samples and views (cell types)
    based on cell number and sample number.

    Args:
        metadata (pd.DataFrame): DataFrame containing sample and view information.
        min_cell_number_per_sample (int): Minimum number of cells required per sample.
        min_sample_per_view (int): Minimum number of samples required per view.
        sample_column (str): Column name containing sample information.
        view_column (str): Column name containing view (cell type) information.

    Returns:
        tuple: A tuple containing two elements. The first element is an array of sample views,
        and the second element is a dictionary containing available samples for each view.

    """

    if (metadata[sample_column].dtype == 'category') & (metadata[view_column].dtype == 'category'):
        dataframe1 = pd.crosstab(metadata[sample_column],metadata[view_column])
        dataframe1 = dataframe1 >= min_cell_number_per_sample
        result_dict = {}
        #filter_sample_below_min_cell
        for column in dataframe1.columns:
            samples = dataframe1.index[dataframe1[column]].tolist()
            result_dict[column] = samples
        #filter_cells_without_enough_samples
        result_dict = {cell_type: samples for cell_type, samples in result_dict.items() if len(samples) >= min_sample_per_view}
        output1 = np.array(list(result_dict.keys()))
        output2 = result_dict
        print('Done_Qc')
    else:
        raise ValueError("Error: The columns you provided are not both category data types")
    return output1,output2

#calulate and generate mean or median of PAS each sample*view

def pseudobulk_matrix(
    metadata:pd.DataFrame = None,
    PAS_dataframe:pd.DataFrame = None,
    usedict:dict = None,
    sample_column:str = None,
    view_column:str = None):
    """
    Generate sample*view pseudobulk matrices based on the given metadata and PAS data.

    Args:
        metadata (pd.DataFrame): DataFrame containing sample and view information.
        PAS_dataframe (pd.DataFrame): DataFrame containing PAS data.
        usedict (dict): Dictionary containing sample view information.
        sample_column (str): Column name containing sample information.
        view_column (str): Column name containing view (cell type) information.

    Returns:
        dict: A dictionary containing sample*view matrices for each view.
    """

    usemetadata = pd.DataFrame({sample_column:metadata[sample_column].astype(str),view_column:metadata[view_column].astype(str)})
    for key in usedict['view_sample'].keys():
        usemetadata_view = usemetadata[usemetadata[view_column] == key].copy()
        usemetadata_view = usemetadata_view[usemetadata_view[sample_column].isin(usedict['view_sample'][key])].copy()
        usedataframe2 = PAS_dataframe.loc[usemetadata_view.index].copy()
        usedataframe2[sample_column] = usemetadata_view[sample_column].copy()
        view_dataframe = usedataframe2.groupby(sample_column,observed=False).mean()
        variance_per_column = view_dataframe.var()
        non_zero_variance_columns = view_dataframe.columns[variance_per_column != 0]
        view_dataframe = view_dataframe[non_zero_variance_columns]
        usedict['view_sample_matrix'][key] = view_dataframe.copy()
    print('Done_matrix')
    return usedict

#generate long table for each matrix because mofa require such input

def change_to_long(usedict:dict = None):
   
    """
    Convert sample*view matrices to long-format tables to meet MOFA's input requirements.

    Args:
        usedict (dict): Dictionary containing sample view information.

    Returns:
        dict: A dictionary containing the converted long-format tables for each view.
    """

    for key in usedict['view_sample_matrix'].keys():
        view_dataframe_long = usedict['view_sample_matrix'][key].copy()
        view_dataframe_long['sample'] = view_dataframe_long.index
        view_dataframe_long = view_dataframe_long.melt(id_vars='sample')
        view_dataframe_long['view'] = key
        view_dataframe_long.rename(columns={'variable': 'feature'}, inplace=True)
        usedict['view_sample_long'][key] = view_dataframe_long
    print('Done_longtable')
    return usedict

# generate a dict countaining all information

def generate_scPAFA_dict(
    metadata:pd.DataFrame = None,
    PAS_dataframe:pd.DataFrame = None,
    min_cell_number_per_sample: int = 10,
    min_sample_per_view: int = 15,
    sample_column:str = None,
    view_column:str = None):

    """
    Generate a dictionary containing all information for scPAFA, including PAS data, sample*view matrices,
    and long-format tables.

    Args:
        metadata (pd.DataFrame): adata.obs,dataFrame containing sample and view information.
        PAS_dataframe (pd.DataFrame): row as cell,column as pathway,dataFrame containing PAS data(the output of pyUCell or other method).
        min_cell_number_per_sample (int): Minimum number of cells required per sample.
        min_sample_per_view (int): Minimum number of samples required per view.
        sample_column (str): Column name containing sample information.
        view_column (str): Column name containing view (cell type) information.

    Returns:
        dict: A dictionary containing PAS data, sample*view matrices, and long-format tables.
    """

    if metadata is None or not isinstance(metadata, pd.DataFrame):
        raise ValueError("Error: 'metadata' must be a valid DataFrame")
    
    if PAS_dataframe is None or not isinstance(PAS_dataframe, pd.DataFrame):
        raise ValueError("Error: 'PAS_dataframe' must be a valid DataFrame")
        
    if sample_column not in metadata.columns or view_column not in metadata.columns:
        raise ValueError("Error: 'sample_column' or 'view_column' not found in metadata")
        
    scPAFA_dict = {}
    scPAFA_dict['PAS'] = PAS_dataframe.copy()
    scPAFA_dict['view'],scPAFA_dict['view_sample'] = pseudobulk_qc(metadata = metadata,
                                                                   sample_column = sample_column,
                                                                   view_column = view_column,
                                                                   min_cell_number_per_sample = min_cell_number_per_sample,
                                                                   min_sample_per_view = min_sample_per_view)
    scPAFA_dict['view_sample_matrix'] ={}
    scPAFA_dict['view_sample_long'] ={}
    scPAFA_dict = pseudobulk_matrix(metadata=metadata,PAS_dataframe=PAS_dataframe,usedict = scPAFA_dict,
                                   sample_column=sample_column,view_column=view_column)
    scPAFA_dict = change_to_long(scPAFA_dict)
    
    long_dataframes_list = [df for df in scPAFA_dict['view_sample_long'].values()]
    combined_long_df = pd.concat(long_dataframes_list, axis=0)
    combined_long_df.reset_index(drop=True, inplace=True)
    scPAFA_dict['long_table_for_mofa'] = combined_long_df
    scPAFA_dict['long_table_for_mofa']['value'] = scPAFA_dict['long_table_for_mofa']['value']*10
    
    return scPAFA_dict
  
