import pandas as pd
from scPAFA.run.entry_point import entry_point

def run_mofapy2(
 long_table:pd.DataFrame,
 scale_groups:bool = False,
 scale_views:bool = False,
 center_groups:bool=True,
 factor_number:int = 10,
 random_seed:int = 0,
 spikeslab_factors:bool = False,
 spikeslab_weights:bool = True,
 ard_factors:bool = False,
 ard_weights:bool = True,
 dropR2:float = None,
 convergence_mode:str ='fast'):
   """
   A wrap function to run the MOFAPY2 and return the MOFA object.
   For detailed parameters and usage, please refer to the MOFAPY2 documentation.
   https://github.com/bioFAM/mofapy2/blob/master/mofapy2/notebooks/getting_started_python.ipynb
   
   Parameters
   ----------
   long_table : pd.DataFrame
      The long-format table containing the data. 
   scale_groups : bool, optional
      Whether to scale the groups.if groups have different ranges/variances, it is good practice to scale each group to unit variance. Default is False.
   scale_views : bool, optional
      Whether to scale the views.if views have different ranges/variances, it is good practice to scale each view to unit variance. Default is False.
   factor_number : int
      The number of factors in the model. Default is 10.
   spikeslab_factors : bool
      Spike-and-slab prior on the factors Default is False.
   spikeslab_weights : bool, optional
      use spike-slab sparsity prior in the weights? Default is TRUE.
   ard_factors :  bool, optional
      use ARD prior in the weights? Default is False. If using multigroup framework, should set to True.
   ard_weights :  bool, optional
      use ARD prior in the weights? Default is TRUE if using multiple views.
   random_seed : int, optional
      Random seed for controlling randomness. Default is 0.
   dropR2 : float, optional
      R2 threshold for factor drop. Minimum variance explained criteria to drop factors while training. Default is None.
   convergence_mode : str, optional
      Convergence mode, can be "fast" (default), "medium", "slow".Default is "fast".

   Returns
   -------
   entry_point
      MOFA object.
   """
   ent = entry_point()
   ent.set_data_options(scale_groups=scale_groups,scale_views = scale_views,center_groups=center_groups)
   ent.set_data_df(long_table.copy()) # Use "gaussian" (default) likelihoods , which are suitable for pseudobulk PAS data
   ent.set_model_options(factors = factor_number,spikeslab_factors=spikeslab_factors,spikeslab_weights = spikeslab_weights,ard_factors = ard_factors,ard_weights =ard_weights)
   ent.set_train_options(convergence_mode = convergence_mode,dropR2 = dropR2,seed = random_seed)
   ent.build()
   ent.run()
   return ent
