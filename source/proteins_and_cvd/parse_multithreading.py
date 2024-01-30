# %%
import os
import pandas as pd
from pathlib import Path


path = "results/incremental_parallel/hosp/agesex_interaction/"
csv_folder = Path(path)  # path to your folder, e.g. to `2022`

df = pd.concat(pd.read_csv(p) for p in csv_folder.glob('**/*.csv'))
df = df.iloc[:, 1:]
df.to_csv(path + "/merged_results.csv")

bonf_threshold = 0.05 / 439
significant = df.loc[df['p'] < bonf_threshold]
# remember about cox assumptions! I didnt check them here
significant.to_csv(path + "/merged_results_bonf_significant.csv")

full_and_significant = significant[significant["formula"] == ("age+sex*protein+avg_sys+"
                                                              "Total_cholesterol+HDL_cholesterol+pack_years"
                                                              "+rheum_arthritis_Y+diabetes_Y+years+rank+on_pill")]
full_and_significant.to_csv(path + "/merged_results_bonf_significant_full.csv")

events = significant["event"].unique()

for event in events:
    event_spec_df = significant.loc[significant["event"] == event]
    event_spec_df.to_csv(path + f"/{event}/{event}_significant.csv")
