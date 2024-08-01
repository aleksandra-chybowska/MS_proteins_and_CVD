# %%
import os
import pandas as pd
from pathlib import Path


path = "results/incremental_parallel/hosp/agesex_interaction/40-69/"
csv_folder = Path(path)

df = pd.concat(pd.read_csv(p) for p in csv_folder.glob('*/*.csv'))
df = df.iloc[:, 1:]
df.to_csv(path + "/merged_results.csv", index=False)

#%%
bonf_threshold = 0.05 / 439
significant = df.loc[df['p'] < bonf_threshold]
# remember about cox assumptions! I didnt check them here
significant.to_csv(path + "/merged_results_bonf_significant.csv", index=False)

# full_formula = ("age+sex*protein+avg_sys+Total_cholesterol+HDL_cholesterol+pack_years+"
#                 "rheum_arthritis_Y+diabetes_Y+years+rank+on_pill")
full_formula = ("age+sex*protein+avg_sys+Total_cholesterol+HDL_cholesterol+pack_years+"
                "rheum_arthritis_Y+diabetes_Y+years+rank+on_pill")
full_and_significant = significant[significant["formula"] == full_formula]
full_and_significant.to_csv(path + "/merged_results_bonf_significant_full.csv")
