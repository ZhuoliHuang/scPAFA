import pandas as pd
import mofax as mfx
from scipy import stats
from scipy.stats import mannwhitneyu, wilcoxon, kruskal
from scipy.stats import ttest_rel, ttest_ind, f_oneway
from scipy.stats import pearsonr, spearmanr, kendalltau
from statsmodels.stats.multitest import fdrcorrection

def normality_and_variance_homogeneity(
        sample_factor_df:pd.DataFrame = None,
        metadata:pd.DataFrame = None, 
        label_column:str = None):
    """
    Perform normality and variance homogeneity tests for each factor and group.

    Parameters:
    -----------
    sample_factor_df : pd.DataFrame
        DataFrame containing sample factors, where rows are samples, columns are factors, and values are floats. Generated from mofa model.
    metadata : pd.DataFrame
        DataFrame containing metadata including the grouping information.The index of the dataframe must contain all samples in mofa model.
    label_column : str
        Name of the column in 'metadata' to use for grouping.

    Returns:
    --------
    result_df : pd.DataFrame
        DataFrame containing test results for each factor and group.
    """    
    # Check if the indices of the dataframes are equal
    if not sample_factor_df.index.equals(metadata.index):
        raise ValueError("Indices of the sample_factor_df and metadata are not equal")

    # Create a result DataFrame to store the test results
    result_df = pd.DataFrame(columns=['Factor', 'Group', 'Normality', 'Normality Test Type', 'Variance Homogeneity', 'Test Type'])

    # Get all the factor columns
    factors = sample_factor_df.columns

    for factor in factors:
        # Get the values for the current factor
        factor_values = sample_factor_df[factor]

        # Get unique groups
        groups = metadata[label_column].unique()

        for group in groups:
            # Filter samples based on the group
            group_samples = factor_values[metadata[label_column] == group]

            # Normality test (Shapiro-Wilk test)
            _, normality_p_value = stats.shapiro(group_samples)

            # Variance homogeneity test (Levene test)
            other_groups = [g for g in groups if g != group]
            other_group_samples = [factor_values[metadata[label_column] == g] for g in other_groups]

            _, variance_homogeneity_p_value = stats.levene(*other_group_samples, group_samples)
            test_type = 'Levene'

            # Create a DataFrame for each result
            each_result = pd.DataFrame({'Factor': factor, 'Group': group, 'Normality': normality_p_value,
                                         'Normality Test Type': 'Shapiro-Wilk',
                                         'Variance Homogeneity': variance_homogeneity_p_value, 'Test Type': test_type}, index=[0])
            
            # Append the result to the result DataFrame
            result_df = pd.concat([result_df, each_result], ignore_index=True)

    # Check if both normality and variance homogeneity tests pass
    if all(result_df['Normality'].values > 0.05) and all(result_df['Variance Homogeneity'].values > 0.05):
        print('Normality and Variance Homogeneity test passed, parametric_test_category() is recommended' )
    else:
        print('Normality or Variance Homogeneity test failed, nonparametric_test_category() is recommended')

    return result_df

def nonparametric_test_category(sample_factor_df, metadata, label_column, pair=False):
    """
    Calculate the significance of differences between different groups for each factor using nonparametric tests 
    (Mann-Whitney U test, Wilcoxon signed-rank test, or Kruskal-Wallis test) and perform Benjamini-Hochberg correction.

    Parameters:
    -----------
    sample_factor_df : pd.DataFrame
        DataFrame containing sample factors, where rows are samples, columns are factors, and values are floats. Generated from mofa model.
    metadata : pd.DataFrame
        DataFrame containing metadata including the grouping information.The index of the dataframe must contain all samples in mofa model.
    label_column : str
        Name of the column in 'metadata' to use for grouping.
    pair : bool, optional
        If True, use Wilcoxon signed-rank test for two-group comparisons; otherwise, use Mann-Whitney U test
        (default is False).

    Returns:
    --------
    result_df : pd.DataFrame
        DataFrame containing the P-values and Benjamini-Hochberg adjusted P-values for each factor.
    """
    # Check if the indices of the dataframes are equal
    if not sample_factor_df.index.equals(metadata.index):
        raise ValueError("Indices of the sample_factor_df and metadata are not equal")
    
    # Create an empty result DataFrame
    result_df = pd.DataFrame(columns=['Factor', 'p_value'])

    # Iterate through each column of sample factors
    for factor in sample_factor_df.columns:
        # Get the values of the current factor
        factor_values = sample_factor_df[factor]

        # Get the classification data corresponding to the factor
        categories = metadata[label_column]

        # Check the number of unique values in the classification
        unique_categories = categories.unique()

        # Choose the appropriate nonparametric test based on the number of categories
        if len(unique_categories) == 2:
            # Only two groups
            if pair:
                # Use the Wilcoxon test
                group1 = factor_values[categories == unique_categories[0]]
                group2 = factor_values[categories == unique_categories[1]]
                _, p_value = wilcoxon(group1, group2)
                test_method = 'Wilcoxon signed-rank test'
            else:
                # Use the Mann-Whitney U test
                group1 = factor_values[categories == unique_categories[0]]
                group2 = factor_values[categories == unique_categories[1]]
                _, p_value = mannwhitneyu(group1, group2)
                test_method = 'Mann-Whitney U test'
        else:
            # Multiple groups, use the Kruskal-Wallis test
            groups = [factor_values[categories == category] for category in unique_categories]
            _, p_value = kruskal(*groups)
            test_method = 'Kruskal-Wallis test'
        
        # Add the result to the result DataFrame
        each_result = pd.DataFrame({'Factor': factor, 'p_value': p_value}, index=[0])
        result_df = pd.concat([result_df, each_result], ignore_index=True)
    
    # Perform Benjamini-Hochberg correction
    _, bh_adjusted_p_values = fdrcorrection(result_df['p_value'])
    result_df['p_adj'] = bh_adjusted_p_values
    result_df['method'] = test_method
    return result_df


def parametric_test_category(sample_factor_df, metadata, label_column, pair=False):
    """
    Calculate the significance of differences between different groups for each factor using parametric tests 
    (T-test, paired T-test, or one-way ANOVA) and perform Benjamini-Hochberg correction.

    Parameters:
    -----------
    sample_factor_df : pd.DataFrame
        DataFrame containing sample factors, where rows are samples, columns are factors, 
        and values are floats.
    metadata : pd.DataFrame
        DataFrame containing sample classifications.
    label_column : str
        The name of the label column in the metadata DataFrame.
    pair : bool, optional
        If True, use paired T-test for two-group comparisons; otherwise, use T-test
        (default is False).

    Returns:
    --------
    result_df : pd.DataFrame
        DataFrame containing the P-values and Benjamini-Hochberg adjusted P-values for each factor.
    """
    
    # Check if the indices of the dataframes are equal
    if not sample_factor_df.index.equals(metadata.index):
        raise ValueError("Indices of the sample_factor_df and metadata are not equal")

    # Create an empty result DataFrame
    result_df = pd.DataFrame(columns=['Factor', 'p_value'])

    # Iterate through each column of sample factors
    for factor in sample_factor_df.columns:
        # Get the values of the current factor
        factor_values = sample_factor_df[factor]

        # Get the classification data corresponding to the factor
        categories = metadata[label_column]

        # Check the number of unique values in the classification
        unique_categories = categories.unique()

        # Choose the appropriate parametric test based on the number of categories
        if len(unique_categories) == 2:
            # Only two groups
            if pair:
                # Use paired T-test
                group1 = factor_values[categories == unique_categories[0]]
                group2 = factor_values[categories == unique_categories[1]]
                _, p_value = ttest_rel(group1, group2)
                test_method = 'paired T-test'
            else:
                # Use T-test
                group1 = factor_values[categories == unique_categories[0]]
                group2 = factor_values[categories == unique_categories[1]]
                _, p_value = ttest_ind(group1, group2)
                test_method = 'T-test'
        else:
            # Multiple groups, use one-way ANOVA
            groups = [factor_values[categories == category] for category in unique_categories]
            _, p_value = f_oneway(*groups)
            test_method = 'one-way ANOVA'
        
        # Add the result to the result DataFrame
        each_result = pd.DataFrame({'Factor': factor, 'p_value': p_value}, index=[0])
        result_df = pd.concat([result_df, each_result], ignore_index=True)
    
    # Perform Benjamini-Hochberg correction
    _, bh_adjusted_p_values = fdrcorrection(result_df['p_value'])
    result_df['p_adj'] = bh_adjusted_p_values
    result_df['method'] = test_method
    return result_df

import pandas as pd
from scipy.stats import pearsonr, spearmanr, kendalltau

def cal_correlation(sample_factor_df, query_df,
                          method:str = 'pearson'):
    
    """
    Calculate correlation coefficients and p-values between pairs of columns from two pandas DataFrames.

    Parameters:
        sample_factor_df : pd.DataFrame
            DataFrame containing sample factors, where rows are samples, columns are factors, 
            and values are floats.
        query_df :pd.DataFrame
            The second DataFrame for which correlations will be calculated.
        method : str
            The method used for calculating correlation. Default is 'pearson', but it can
            also be 'spearman' or 'kendall'.

    Returns:
        DataFrame: A DataFrame containing the correlation coefficients and p-values for each pair of columns.
            Columns in the result DataFrame include 'column1', 'column2', 'Correlation','p_value' and 'method'.

    Raises:
        ValueError: If the indices of the two DataFrames are not equal or if an unsupported correlation method is specified.
    """
    
    # Check if the indices of the dataframes are equal
    if not sample_factor_df.index.equals(query_df.index):
        raise ValueError("Indices of the dataframes are not equal")

    # Initialize an empty result dataframe
    result_df = pd.DataFrame(columns=['column1', 'column2', 'Correlation', 'P-Value'])

    # Get the column names of dataframes A and B
    a_columns = sample_factor_df.columns
    b_columns = query_df.columns

    # Iterate through each pair of columns in dataframes A and B
    for a_col in a_columns:
        for b_col in b_columns:
            # Extract the data from the columns
            a_data = sample_factor_df[a_col]
            b_data = query_df[b_col]

            # Calculate the correlation and p-value based on the selected method
            if method == 'pearson':
                correlation, p_value = pearsonr(a_data, b_data)
            elif method == 'spearman':
                correlation, p_value = spearmanr(a_data, b_data)
            elif method == 'kendall':
                correlation, p_value = kendalltau(a_data, b_data)
            else:
                raise ValueError("Unsupported correlation method")

            # Add the results to the result dataframe
            each_result = pd.DataFrame({'column1': a_col, 'column2': b_col, 'Correlation': correlation, 'P-Value': p_value}, index=[0])
            result_df = pd.concat([result_df, each_result], ignore_index=True)
    
    result_df['method'] = method
    return result_df