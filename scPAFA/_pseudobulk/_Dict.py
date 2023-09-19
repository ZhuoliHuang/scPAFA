import numpy as np
import pandas as pd

# A QC funtion for selecting usable sample and view(celltype),based on cell number and sample number

def Pseudobulk_Qc(
    metadata:pd.DataFrame = None,
    min_cell_number_per_sample: int = 10,
    min_sample_per_view: int = 15,
    sample_column:str = None,
    view_column:str = None):
    if (metadata[sample_column].dtype == 'category') & (metadata[view_column].dtype == 'category'):
        dataframe1 = pd.crosstab(metadata[sample_column],metadata[view_column])
        dataframe1 = dataframe1 >= min_cell_number_per_sample
        result_dict = {}
        for column in dataframe1.columns:
            samples = dataframe1.index[dataframe1[column]].tolist()
            result_dict[column] = samples
            result_dict = {cell_type: samples for cell_type, samples in result_dict.items() if len(samples) >= min_sample_per_view}
        output1 = np.array(list(result_dict.keys()))
        output2 = result_dict
        print('Done_Qc')
    else:
        raise ValueError("Error: The columns you provided are not both category data types")
    return output1,output2

#calulate and generate mean or median of PAS each sample*view

def Pseudobulk_Matrix(
    metadata:pd.DataFrame = None,
    PAS_dataframe:pd.DataFrame = None,
    usedict:dict = None,
    sample_column:str = None,
    view_column:str = None):
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

def Change_to_Long(usedict:dict = None):
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

def Generate_scPAFA_Dict(
    metadata:pd.DataFrame = None,
    PAS_dataframe:pd.DataFrame = None,
    min_cell_number_per_sample: int = 10,
    min_sample_per_view: int = 15,
    sample_column:str = None,
    view_column:str = None):
    
    if metadata is None or not isinstance(metadata, pd.DataFrame):
        raise ValueError("Error: 'metadata' must be a valid DataFrame")
    
    if PAS_dataframe is None or not isinstance(PAS_dataframe, pd.DataFrame):
        raise ValueError("Error: 'PAS_dataframe' must be a valid DataFrame")
        
    if sample_column not in metadata.columns or view_column not in metadata.columns:
        raise ValueError("Error: 'sample_column' or 'view_column' not found in metadata")
        
    scPAFA_dict = {}
    scPAFA_dict['PAS'] = PAS_dataframe.copy()
    scPAFA_dict['view'],scPAFA_dict['view_sample'] = Pseudobulk_Qc(metadata = metadata,
                                                                   sample_column = sample_column,
                                                                   view_column = view_column,
                                                                   min_cell_number_per_sample = min_cell_number_per_sample,
                                                                   min_sample_per_view = min_sample_per_view)
    scPAFA_dict['view_sample_matrix'] ={}
    scPAFA_dict['view_sample_long'] ={}
    scPAFA_dict = Pseudobulk_Matrix(metadata=metadata,PAS_dataframe=PAS_dataframe,usedict = scPAFA_dict,
                                   sample_column=sample_column,view_column=view_column)
    scPAFA_dict = Change_to_Long(scPAFA_dict)
    
    long_dataframes_list = [df for df in scPAFA_dict['view_sample_long'].values()]
    combined_long_df = pd.concat(long_dataframes_list, axis=0)
    combined_long_df.reset_index(drop=True, inplace=True)
    scPAFA_dict['long_table_for_mofa'] = combined_long_df
    scPAFA_dict['long_table_for_mofa']['value'] = scPAFA_dict['long_table_for_mofa']['value']*10
    
    return scPAFA_dict
  
