import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from math import ceil

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

    """

    if (metadata[sample_column].dtype == 'category') & (metadata[view_column].dtype == 'category'):
        dataframe1 = pd.crosstab(metadata[sample_column],metadata[view_column])
        dataframe1 = dataframe1 >= min_cell_number_per_sample
        result_dict = {}
        #filter_sample_below_min_cell
        for column in dataframe1.columns:
            samples = dataframe1.index[dataframe1[column]].tolist()
            result_dict[column] = samples
        #filter_views_without_enough_samples
        result_dict = {cell_type: samples for cell_type, samples in result_dict.items() if len(samples) >= min_sample_per_view}
        output1 = np.array(list(result_dict.keys()))
        output2 = result_dict
        print('Done_Qc')
    else:
        raise ValueError("Error: The columns you provided are not all category data types")
    return output1,output2

#calulate and generate mean or median of PAS each sample*view

def pseudobulk_matrix(
    metadata:pd.DataFrame = None,
    PAS_dataframe:pd.DataFrame = None,
    usedict:dict = None,
    sample_column:str = None,
    view_column:str = None,
    number_of_features:float = None,
    regress_out:bool = False,
    formula_use:str = None,
    sample_metadata:pd.DataFrame = None):

    """
    Generate sample*view pseudobulk matrices based on the given metadata and PAS data.

    """

    usemetadata = pd.DataFrame({sample_column:metadata[sample_column].astype(str),view_column:metadata[view_column].astype(str)})
    for key in usedict['view_sample'].keys():
        
        #pseudobulk_mean
        usemetadata_view = usemetadata[usemetadata[view_column] == key].copy()
        usemetadata_view = usemetadata_view[usemetadata_view[sample_column].isin(usedict['view_sample'][key])].copy()
        usedataframe2 = PAS_dataframe.loc[usemetadata_view.index].copy()
        usedataframe2[sample_column] = usemetadata_view[sample_column].copy()
        view_dataframe = usedataframe2.groupby(sample_column,observed=False).mean()

        #delete the column that variance == 0
        variance_per_column = view_dataframe.var()
        non_zero_variance_columns = view_dataframe.columns[variance_per_column != 0]
        view_dataframe = view_dataframe[non_zero_variance_columns]

        #regressout
        if regress_out:
            sample_metadata_use = sample_metadata.loc[view_dataframe.index,:].copy()
            for pathway in view_dataframe.columns:
                sample_metadata_use['score'] = view_dataframe[pathway]
                model_GLM = smf.glm(formula=formula_use, data=sample_metadata_use,family=sm.families.Gaussian()).fit()
                view_dataframe[pathway] = model_GLM.resid_response
        
        #top_variance
        variance_per_column = view_dataframe.var()
        top_columns = variance_per_column.nlargest(number_of_features).index
        view_dataframe = view_dataframe[top_columns]
        
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
        #avoid sample feature name in different views
        view_dataframe_long.feature = view_dataframe_long.feature.astype(str)+'_View_'+view_dataframe_long.view.astype(str)
        usedict['view_sample_long'][key] = view_dataframe_long
    print('Done_longtable')
    return usedict

# generate a dict countaining all information for mofapy2
def generate_scpafa_input(
    metadata:pd.DataFrame,
    PAS_dataframe:pd.DataFrame,
    sample_column:str,
    view_column:str,
    min_cell_number_per_sample: int = 10,
    min_sample_per_view: int = 15,
    top_percentage:float = 0.25,
    regress_out:bool = False,
    sample_metadata:pd.DataFrame=None,
    regress_categorical_variable:list = [],
    regress_continuous_variable:list = [],
    return_full_dict:bool = False
    ):
    """
    Generate a dictionary containing all information for scPAFA, including sample*view matrices and long-format tables.

    Parameters
    ----------
    metadata : pd.DataFrame
        A DataFrame containing sample and view information, row index as cell. For example: adata.obs.
    PAS_dataframe : pd.DataFrame
        A DataFrame with cells as rows, pathways as columns, containing PAS data.
        (the output of pyUCell or other method).
    sample_column : str
        Column name containing sample information.
    view_column : str
        Column name containing view (cell type) information.
    min_cell_number_per_sample : int
        Minimum number of cells required per sample. Default is 10.
    min_sample_per_view : int
        Minimum number of samples required per view. Default is 15.
    top_percentage: float
        A float >0 and <= 1. Select the top 25%(Default) pathways with the maximum variance across samples.
        (after regress out if regress out == True)
    regress_out : bool
        Wheather to regress out (mostly) unwanted sources of variation in each view by samples uses simple linear regression(GLM).Default is False.
    sample_metadata : pd.DataFrame,
        A Dataframe containing information of samples (should contain all samples in metadata),
        row index as samples, columns as information (for example : batch).
    regress_categorical_variable: list,
        A list of columns in sample_metadata that needs to be regress_out, dtype as category.
    regress_continuous_variable: list ,
        A list of columns in sample_metadata that needs to be regress_out, dtype as int or float.

    Returns
    -------
    dict
        A dictionary containing PAS data, sample*view matrices, and long-format tables.
    """

    if metadata is None or not isinstance(metadata, pd.DataFrame):
        raise ValueError("Error: 'metadata' must be a valid DataFrame")
    
    if PAS_dataframe is None or not isinstance(PAS_dataframe, pd.DataFrame):
        raise ValueError("Error: 'PAS_dataframe' must be a valid DataFrame")
    
    if not all((metadata.index == PAS_dataframe.index)):
        raise ValueError("Error: The index of metadata and PAS_dataframe must be identical")
        
    if sample_column not in metadata.columns or view_column not in metadata.columns:
        raise ValueError("Error: 'sample_column' or 'view_column' not found in metadata")
    
    
    if (top_percentage < 0) or (top_percentage > 1):
        raise ValueError("Error: 'top_percentage' must > 0 and <= 1")

    # how many pathways should be used in mofapy2
    number_of_features = ceil(PAS_dataframe.shape[1]*top_percentage)
    print('Select the top '+str(number_of_features)+' pathways with the maximum variance')

    #generate_regress_out_formula
    formula_regress = []
    if regress_out:
        if sample_metadata is None:
            raise ValueError("Error: no sample_metadata was provided")
        elif (regress_categorical_variable == []) and (regress_continuous_variable == []):
            raise ValueError("Error: no variables for regress were provided")
        else:
            regress_categorical_variable = ["C(" + element + ")" for element in regress_categorical_variable]
            formula_regress.extend(regress_categorical_variable)
            formula_regress.extend(regress_continuous_variable)
            formula_regress = "+".join(formula_regress)
            formula_regress = 'score~'+formula_regress
            print('The regress formula is ' + formula_regress)

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
                                   sample_column=sample_column,view_column=view_column,number_of_features=number_of_features,
                                   regress_out=regress_out,formula_use=formula_regress,sample_metadata=sample_metadata)
    
    scPAFA_dict = change_to_long(scPAFA_dict)
    
    long_dataframes_list = [df for df in scPAFA_dict['view_sample_long'].values()]
    combined_long_df = pd.concat(long_dataframes_list, axis=0)
    combined_long_df.reset_index(drop=True, inplace=True)
    scPAFA_dict['long_table_for_mofa'] = combined_long_df
    del scPAFA_dict['PAS']
    
    if return_full_dict:
        return scPAFA_dict
    else:
        return scPAFA_dict['long_table_for_mofa']
  
# generrate a multi group pseudobulk dataframe to run mofapy2
def generate_scpafa_input_multigroup(
    metadata:pd.DataFrame,
    PAS_dataframe:pd.DataFrame,
    sample_column:str,
    view_column:str,
    group_column:str,
    min_cell_number_per_sample: int = 10,
    min_percentage_sample_per_view: float = 0.75,
    min_sample_per_view:int = 15,
    top_percentage:float = 0.25):

    """
    Generate a dictionary containing all information for scPAFA, including sample*view matrices and long-format tables.

    Parameters
    ----------
    metadata : pd.DataFrame
        A DataFrame containing sample and view information, row index as cell. For example: adata.obs.
    PAS_dataframe : pd.DataFrame
        A DataFrame with cells as rows, pathways as columns, containing PAS data.
        (the output of pyUCell or other method).
    sample_column : str
        Column name containing sample information.
    view_column : str
        Column name containing view (cell type) information.
    min_cell_number_per_sample : int
        Minimum number of cells required per sample. Default is 10.
    min_percentage_sample_per_view:
        The percentage of minimum number of samples required per view in each group. Default is 0.75. 
    min_sample_per_view : int
        Minimum number of samples required per view. Default is 15. 
    top_percentage: float
        A float >0 and <= 1. Select the top 25%(Default) pathways with the maximum variance across samples.

    Returns
    -------
    dict
        A long-format tables.
    """

    if metadata is None or not isinstance(metadata, pd.DataFrame):
        raise ValueError("Error: 'metadata' must be a valid DataFrame")
    
    if PAS_dataframe is None or not isinstance(PAS_dataframe, pd.DataFrame):
        raise ValueError("Error: 'PAS_dataframe' must be a valid DataFrame")
    
    if not all((metadata.index == PAS_dataframe.index)):
        raise ValueError("Error: The index of metadata and PAS_dataframe must be identical")

    if (min_percentage_sample_per_view < 0) or (min_percentage_sample_per_view > 1):
        raise ValueError("Error: 'min_percentage_sample_per_view' must > 0 and <= 1")

    if sample_column not in metadata.columns or view_column not in metadata.columns or group_column not in metadata.columns:
        raise ValueError("Error: 'sample_column' or 'view_column' or 'group_column' not found in metadata")
    
    if (top_percentage < 0) or (top_percentage > 1):
        raise ValueError("Error: 'top_percentage' must > 0 and <= 1")

    if (metadata[sample_column].dtype == 'category') & (metadata[group_column].dtype == 'category'):
        dataframe1 = pd.crosstab(metadata[sample_column],metadata[group_column])
    else:
        raise ValueError("Error: The columns you provided are not all category data types")
    
    if len(metadata[group_column].value_counts()) < 2:
        raise ValueError("Error: Groups should more than 2")
    else:
        print(str(len(metadata[group_column].value_counts()))+ ' groups indentified')
        
    count_non_zero = dataframe1[dataframe1 != 0].count(axis=1)
    if count_non_zero.sum() > len(count_non_zero):
        raise ValueError("Error: Different groups have overlap samples")
    
    # to_cal_number_of_samples_in_each_group
    sample_meta = metadata[[sample_column,group_column]]
    sample_meta = sample_meta.drop_duplicates(subset=sample_column,ignore_index=True)
    sample_series = sample_meta[group_column].value_counts()*min_percentage_sample_per_view
    
    longdf_list = []
    
    for group in metadata[group_column].value_counts().index:
        print('processing group '+group)
        
        group_metadata = metadata[metadata[group_column] == group]
        group_PAS = PAS_dataframe.loc[group_metadata.index]
        
        group_longdf = generate_scpafa_input(metadata=group_metadata,
                                PAS_dataframe=group_PAS,
                                min_cell_number_per_sample=min_cell_number_per_sample,
                                min_sample_per_view=max(min_sample_per_view,round(sample_series[group])),
                                sample_column=sample_column,
                                view_column=view_column,
                                top_percentage=top_percentage,
                                return_full_dict = False)
        
        group_longdf.loc[:,'group'] = group
        
        longdf_list.append(group_longdf)
    
    #multigroup_dataframe
    result_df = pd.concat(longdf_list, axis=0, ignore_index=True)
    
    return result_df