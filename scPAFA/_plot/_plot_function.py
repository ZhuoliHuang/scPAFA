import mofax as mfx
import umap
import pandas as pd
from plotnine import ggplot,geom_point,aes,ggtitle,theme,element_text,element_rect,element_line,labs,scale_color_manual
import seaborn as sns
import random


def runumap_and_plot(
    mofamodel:mfx.core.mofa_model = None,
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
    show:bool = True,
    color_mapping:dict = None,
    retun_umap_embedding:bool = False):
    
    """
    Run UMAP (Uniform Manifold Approximation and Projection) dimensionality reduction
    and create a visualization plot.

    Parameters:
    -----------
    mofamodel : mofax.core.mofa_model
        MOFA model containing factor data to be used for UMAP.
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
    show : bool, optional
        Whether to display the plot (True) or return the plot object (False).
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
    
    factor_dataframe = mofamodel.get_factors(df=True)
    factor_dataframe = factor_dataframe.loc[metadata.index]
    
    reducer = umap.UMAP(n_neighbors=n_neighbors,min_dist=min_dist,random_state=random_state)
    embedding = reducer.fit_transform(factor_dataframe)
    
    umap_dataframe = pd.DataFrame({'UMAP1':embedding[:,0],"UMAP2":embedding[:,1]})
    umap_dataframe.index = factor_dataframe.index
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
    if show:
        print(p)
    
    if retun_umap_embedding:
        return umap_dataframe
    else:
        return p

#some helper function
def text_map_sig(p_value):
    if p_value < 0.01:
        return '**'  # 深色
    elif p_value < 0.05:
        return '*'  # 深色
    else:
        return 'ns'  # 浅色

def color_map_sig(p_value):
    if p_value < 0.01:
        return '#0C7489'  # 深色
    elif p_value < 0.05:
        return '#119DA4'  # 深色
    else:
        return '#F4FAFF'  # 浅色

def assign_colors_to_categories(dataframe, column_name):
    # 获取数据框中唯一的类别
    categories = dataframe[column_name].unique()
    
    # 准备20种不同的低饱和度颜色
    low_saturation_colors = ['#363636', '#1d4d67', '#365952', '#773a00', '#3d2a50',
                        '#575341', '#4d2b25', '#5a5252', '#6d6d6d', '#7e6f59',
                        '#605f14', '#138083', '#b21c43', '#363636', '#605f14',
                        '#363636', '#3d2331', '#2e6d3a', '#3d2a50', '#6b664a']
    # 随机分配颜色给每个类别（不放回抽样）
    color_dict = {}
    sampled_colors = random.sample(low_saturation_colors, len(categories))
    for i, category in enumerate(categories):
        # 将类别和颜色对应存储在字典中
        color_dict[category] = sampled_colors[i]
    
    return color_dict



def draw_cluster_heatmap_category(
    sample_factor_df:pd.DataFrame = None,
    metadata:pd.DataFrame = None,
    label_column:str = None,
    p_value_dataframe:pd.DataFrame = None,
    *args, **kwargs):
    
    """
    Create a clustered heatmap with metadata coloring and significance legend.

    This function generates a clustered heatmap using Seaborn with custom coloring based on
    significance and metadata. It also adds legends for significance and metadata labels.

    Parameters:
    sample_factor_df : pd.DataFrame
        DataFrame containing sample factors, where rows are samples, columns are factors, 
        and values are floats.
    metadata : pd.DataFrame
        DataFrame containing sample classifications.
    label_column :str
        The name of the column in the metadata DataFrame to use for labeling.
    p_value_dataframe : pd.DataFrame
        The DataFrame containing factors and p-values. Generated by stat functions.
    *args, **kwargs: Additional arguments that can be passed to Seaborn's `clustermap` function.

    Returns:
    - sns.ClusterGrid: The Seaborn ClusterGrid object representing the clustered heatmap.
    """
    
    # Check if the indices of the dataframes are equal
    if not sample_factor_df.index.equals(metadata.index):
        raise ValueError("Indices of the sample_factor_df and metadata are not equal")
    
    factor_color_dataframe = p_value_dataframe[['Factor','p_adj']]
    factor_color_dataframe.index = factor_color_dataframe.Factor
    del factor_color_dataframe['Factor']
    
    sig_col = dict({'**':'#0C7489','*':'#119DA4','ns':'#F4FAFF'})
    color_map = assign_colors_to_categories(metadata,label_column)
    factor_color_dataframe['sig'] = factor_color_dataframe['p_adj'].apply(text_map_sig)
    factor_color_dataframe['Color'] = factor_color_dataframe['p_adj'].apply(color_map_sig)
    g = sns.clustermap(sample_factor_df,col_cluster=False,col_colors=factor_color_dataframe['Color'],
                   yticklabels=False,row_colors = metadata[label_column].map(color_map),
                   cmap = sns.diverging_palette(220, 20, as_cmap=True),*args, **kwargs)
    g.fig.subplots_adjust(right=0.6)
    g.ax_cbar.set_position((0.7, .15, .03, .2))
    
    # Draw the legend bar for the classes 
    for label in factor_color_dataframe['sig'].unique():
        g.ax_row_dendrogram.bar(0,0,color=sig_col[label],
                                label=label, linewidth=0)
    g.ax_row_dendrogram.legend(loc=(0,0),bbox_to_anchor=(1.1, 1.1), ncol=3,title = 'Significance',frameon=False)

    for label in metadata[label_column].unique():
        g.ax_col_dendrogram.bar(0,0,color=color_map[label],
                                label=label, linewidth=0)
    g.ax_col_dendrogram.set_position((0.75, .1, .03, .7))
    g.ax_col_dendrogram.legend(loc='best', ncol=1,title = label_column,frameon=False)
    
    return g