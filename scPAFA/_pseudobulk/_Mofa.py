from scPAFA.run.entry_point import entry_point
import pandas as pd

def run_mofapy2(
 long_table:pd.DataFrame = None,
 scale_views:bool = False,
 factor_number:int = 10,
 seed:int = 42,
 dropR2:float = None,
 convergence_mode:str ='fast'):
    ent = entry_point()
    ent.set_data_options(scale_views = scale_views)
    ent.set_data_df(long_table.copy())
    ent.set_model_options(factors = factor_number)
    ent.set_train_options(convergence_mode = convergence_mode,dropR2 = dropR2,seed = seed)
    ent.build()
    ent.run()
    return ent
