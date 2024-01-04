import umap
import pandas as pd
from plotnine import ggplot,geom_point,aes,ggtitle,theme,element_text,element_rect,element_line,labs,scale_color_manual
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from nheatmap import nhm

def runumap_and_plot(
    sample_factor_df:pd.DataFrame = None,
    metadata:pd.DataFrame = None,
    label_column:str = None,
    n_neighbors:int = 15,
    min_dist:float = 0.1,
    random_state:int = 42,
    width:float = 6,
    height:float =5,
    point_size:float = 2,
    title:str = None,
    title_text_size:float = 15,
    axis_text_size:float = 10,
    legend_text_size:float = 10,
    legend_position:str = 'right',
    color_mapping:dict = None,
    retun_umap_embedding:bool = False):
    
    """
    Run UMAP (Uniform Manifold Approximation and Projection) dimensionality reduction
    and create a visualization plot.

    Parameters:
    -----------
    sample_factor_df : pd.DataFrame
        DataFrame containing sample factors, where rows are samples, columns are factors, and values are floats. Generated from mofa model.
    metadata : pd.DataFrame
        Metadata dataframe containing additional information for each data point. The index of the dataframe must contain all samples in mofamodel. 
    label_column : str
        The column in 'metadata' to be used for coloring data points in the plot.
    n_neighbors : int, optional
        The number of neighbors to consider for UMAP.
    min_dist : float, optional
        The minimum distance between points in the UMAP embedding.
    random_state : int, optional
        Random state for reproducibility of UMAP.default = 42
    width : float, optional
        Width of the plot.
    height : float, optional
        Height of the plot.
    point_size : float, optional
        Size of data points in the plot.
    title : str, optional
        Title of the UMAP visualization plot.
    title_text_size : float, optional
        Text size for the plot title.
    axis_text_size : float, optional
        Text size for axis labels.
    legend_text_size : float, optional
        Text size for legend labels.
    legend_position : str, optional
        Position of the legend ('right', 'left', 'top', 'bottom').
    color_mapping:dict, optional
        a dict with key as category and value as color
    retun_umap_embedding:bool ,optional
        decide what ro return 
        default as false,return ggplot object
        if true,return umap embedding dataframe
    Returns:
    --------
    returns the ggplot object representing the UMAP plot or umap embedding pandas dataframe of the UMAP plot.
    """
    sample_factor_df =  sample_factor_df.loc[metadata.index]
    
    reducer = umap.UMAP(n_neighbors=n_neighbors,min_dist=min_dist,random_state=random_state)
    embedding = reducer.fit_transform( sample_factor_df)
    
    umap_dataframe = pd.DataFrame({'UMAP1':embedding[:,0],"UMAP2":embedding[:,1]})
    umap_dataframe.index =  sample_factor_df.index
    umap_dataframe[label_column] = metadata[label_column]
    
    p =(ggplot(aes(x='UMAP1', y='UMAP2',color=label_column), umap_dataframe)
    + geom_point(size=point_size,alpha=1)
    +theme(
        figure_size=[width,height],
        axis_text=element_text(size=axis_text_size,color="black"),
        axis_title_x=element_text(size=axis_text_size,color="black"),
        axis_title_y=element_text(size=axis_text_size,color="black"),
        plot_title=element_text(margin={'b': 1, 'r': 0, 'units': 'pt'},size=title_text_size,color="black",hjust=0.5),
        panel_background=element_rect(fill="white", alpha=0),
        panel_grid_major=element_line(size=0.3, alpha=0.0,color="black"),
        panel_grid_minor=element_line(size=0.3, alpha=0.0,color="black"),
        panel_border=element_rect(color="black", size=1),
        legend_title = element_text(size=legend_text_size),
        legend_text = element_text(size=legend_text_size),
        legend_background=element_rect(size=0.5,alpha=0),
        legend_position=legend_position, 
        legend_direction='vertical',
        legend_key_size=4)
   )
    if color_mapping is not None:
        p = p+scale_color_manual(values=color_mapping)
    if title is not None:
        p = p+ggtitle(title)

    if retun_umap_embedding:
        return umap_dataframe
    else:
        return p
    
def plot_factor_scatter_2D(
    sample_factor_df:pd.DataFrame = None,
    metadata:pd.DataFrame = None,
    label_column:str = None,
    factor_x:str = None,
    factor_y:str = None,
    width:float = 6,
    height:float =5,
    point_size:float = 2,
    title:str = None,
    title_text_size:float = 15,
    axis_text_size:float = 10,
    legend_text_size:float = 10,
    legend_position:str = 'right',
    color_mapping:dict = None):
    
    """
    Run UMAP (Uniform Manifold Approximation and Projection) dimensionality reduction
    and create a visualization plot.

    Parameters:
    -----------
    sample_factor_df : pd.DataFrame
        DataFrame containing sample factors, where rows are samples, columns are factors, and values are floats. Generated from mofa model.
    metadata : pd.DataFrame
        Metadata dataframe containing additional information for each data point. The index of the dataframe must contain all samples in mofamodel. 
    label_column : str
        The column in 'metadata' to be used for coloring data points in the plot.
    factor_x : str
        The factor to be x axis. For example: "Factor1"
    factor_y : str
        The factor to be y axis 
    width : float, optional
        Width of the plot.
    height : float, optional
        Height of the plot.
    point_size : float, optional
        Size of data points in the plot.
    title : str, optional
        Title of the UMAP visualization plot.
    title_text_size : float, optional
        Text size for the plot title.
    axis_text_size : float, optional
        Text size for axis labels.
    legend_text_size : float, optional
        Text size for legend labels.
    legend_position : str, optional
        Position of the legend ('right', 'left', 'top', 'bottom').
    color_mapping:dict, optional
        a dict with key as category and value as color
    retun_umap_embedding:bool ,optional
        decide what ro return 
        default as false,return ggplot object
        if true,return umap embedding dataframe
    Returns:
    --------
    returns the ggplot object representing the UMAP plot or umap embedding pandas dataframe of the UMAP plot.
    """
    sample_factor_df =  sample_factor_df.loc[metadata.index]
    umap_dataframe = pd.DataFrame({factor_x:sample_factor_df.loc[:,factor_x],factor_y:sample_factor_df.loc[:,factor_y]})
    umap_dataframe.index =  sample_factor_df.index
    umap_dataframe[label_column] = metadata[label_column]
    
    p =(ggplot(aes(x=factor_x, y=factor_y,color=label_column), umap_dataframe)
    + geom_point(size=point_size,alpha=1)
    +theme(
        figure_size=[width,height],
        axis_text=element_text(size=axis_text_size,color="black"),
        axis_title_x=element_text(size=axis_text_size,color="black"),
        axis_title_y=element_text(size=axis_text_size,color="black"),
        plot_title=element_text(margin={'b': 1, 'r': 0, 'units': 'pt'},size=title_text_size,color="black",hjust=0.5),
        panel_background=element_rect(fill="white", alpha=0),
        panel_grid_major=element_line(size=0.3, alpha=0.0,color="black"),
        panel_grid_minor=element_line(size=0.3, alpha=0.0,color="black"),
        panel_border=element_rect(color="black", size=1),
        legend_title = element_text(size=legend_text_size),
        legend_text = element_text(size=legend_text_size),
        legend_background=element_rect(size=0.5,alpha=0),
        legend_position=legend_position, 
        legend_direction='vertical',
        legend_key_size=4)
   )
    if color_mapping is not None:
        p = p+scale_color_manual(values=color_mapping)
    if title is not None:
        p = p+ggtitle(title)
    return p

#some helper function
def text_map_sig(p_value):
    if p_value < 0.01:
        return '**'  # 深色
    elif p_value < 0.05:
        return '*'  # 深色
    else:
        return 'ns'  # 浅色

def draw_cluster_heatmap(
    sample_factor_df:pd.DataFrame = None,
    sample_annotaion_df:pd.DataFrame = None,
    p_value_dataframe:pd.DataFrame = None,
    cmapCenter:str='RdBu_r',
    custom_cmap=None,
    *args, **kwargs):
    
    """
    Create a clustered heatmap with metadata coloring and significance legend.

    This function generates a clustered heatmap with custom coloring based on
    significance and metadata. It also adds legends for significance and metadata labels.

    Parameters:
    sample_factor_df : pd.DataFrame
        DataFrame containing sample factors, where rows are samples, columns are factors, 
        and values are floats.
    sample_annotaion_df :
        A pandas dataframe, index as samples,columns as annotation
    p_value_dataframe : pd.DataFrame
        The DataFrame containing factors and p-values. Generated by stat functions.
    *args, **kwargs: Additional arguments that can be passed to nheatmap's `nhm` function.

    Returns:
    - sns.ClusterGrid: The Seaborn ClusterGrid object representing the clustered heatmap.
    """
    
    # Check if the indices of the dataframes are equal
    if not sample_factor_df.index.equals(sample_annotaion_df.index):
        raise ValueError("Indices of the sample_factor_df and sample_annotaion_df are not equal")
    
    factor_color_dataframe = p_value_dataframe[['Factor','p_adj']]
    factor_color_dataframe.index = factor_color_dataframe.Factor
    del factor_color_dataframe['Factor']
    factor_color_dataframe['Significance'] = factor_color_dataframe['p_adj'].apply(text_map_sig)
    del factor_color_dataframe['p_adj']
    factor_color_dataframe=factor_color_dataframe.loc[sample_factor_df.columns]
    
    colors_dict = {'**': '#0C7489', '*': '#119DA4', 'ns': '#F4FAFF'}
    bounds = sorted(set(factor_color_dataframe['Significance'].values))
    cmap_sig = ListedColormap([colors_dict[key] for key in bounds])
    
    dict_sig = {'Significance':cmap_sig}
    if custom_cmap is not None:
        dict_sig = dict_sig | custom_cmap

    g = nhm(data=sample_factor_df,dfr=sample_annotaion_df,dfc=factor_color_dataframe,linewidths=0, cmaps=dict_sig, showyticks=False,cmapCenter=cmapCenter,*args, **kwargs)
    g.hcluster(optimal_ordering=True,col_cluster=False)
    fig, plots = g.run()
    
    return fig

def plot_weights_butterfly(weight_matrix, factor_name, figsize=(10, 8), positive_color='salmon', negative_color='lightblue',
                           n_largest:int = 10,n_smallest:int = 10,label_beside_bar = True):
    
    """
    Create a butterfly plot displaying the top positive and negative weights associated with a given factor.

    Parameters:
    - weight_matrix (pd.DataFrame): DataFrame containing weight values associated with factors. Can be generated by mofax get_weights function. For example,
    weight_matrix = pd.concat(MOFA_model.get_weights(concatenate_views =False,scale=True,df =True))
    - factor_name (str): Name of the factor for which weights are to be plotted.
    - figsize (tuple, optional): Size of the plot figure (width, height). Defaults to (10, 8).
    - positive_color (str, optional): Color for positive weights. Defaults to 'salmon'.
    - negative_color (str, optional): Color for negative weights. Defaults to 'lightblue'.
    - n_largest (int, optional): Top n features to plot.
    - n_smallest (int, optional):  Smallest n features to plot.
    - label_beside_bar(bool, optional): True or False
    Returns:
    - plt: Matplotlib plot object displaying the butterfly plot of top positive and negative weights.
    """

    factor_selected = weight_matrix[factor_name]

    top_10_positive = factor_selected.nlargest(n_largest)
    top_10_negative = factor_selected.nsmallest(n_smallest)

    top_10 = pd.concat([top_10_positive, top_10_negative])
    top_10 = top_10.sort_values(ascending=False)

    plt.figure(figsize=figsize)
    ax = sns.barplot(x=top_10.values, y=top_10.index, palette=[positive_color if x >= 0 else negative_color for x in top_10.values])

    if label_beside_bar:
        for i, v in enumerate(top_10.values):
            if v >= 0:
                ax.text(-0.1, i, str(top_10.index[i]), va='center', ha='right')
            else:
                ax.text(0.1, i, str(top_10.index[i]), va='center')
        ax.yaxis.set_visible(False)
        sns.despine(left=True)
        
    plt.xlabel('Weights')
    plt.ylabel('')
    plt.title(factor_name)
    
    return plt