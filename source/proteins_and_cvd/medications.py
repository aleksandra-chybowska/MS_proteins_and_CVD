import os
import pandas as pd
from lifelines import CoxPHFitter

from lib.cox import summary_and_test
from lib.stats import summary

#%%
flag = "hosp"  # hosp_gp, hosp, hosp_gp_cons
run = "agesex"

meds = pd.read_table('data/disease/medications/all_CVD_prescriptions.txt')
treated = meds["id"].unique()

# filter proteins to significant associations
proteins = pd.read_csv('results/cox/hosp/prepped/proteins_hosp_all_events_scaled_8660.csv')
proteins.set_index("id", inplace=True)
sig = pd.read_csv("results/incremental_parallel/hosp/agesex/merged_results_bonf_significant_full.csv")
proteins = proteins.loc[:, sig["id"]]
annots = pd.read_csv("data/annotations/short_annots.csv")

path = f'results/medications/{run}/{flag}'
if not os.path.exists(path):
    os.makedirs(path)
    print(f"Path: {path} created!")

# look only at events for which significant associations were repored
events = sig["event"].unique()
# %%
for event in events:
    print(event)
    cox_path = f"results/cox/{flag}/prepped/cox_{flag}_{event}_prepped.csv"
    cox = pd.read_csv(cox_path)
    cox["treated"] = 0
    cox.loc[cox["id"].isin(treated), "treated"] = 1
    cox.set_index("id", inplace=True)
    # 1786 treated of n = 8660
    df = pd.merge(cox, proteins, how="inner", left_index=True, right_index=True)

    full = []
    for protein in proteins.columns:
        print(protein)
        cph = CoxPHFitter()
        formula = (f"age+sex+{protein}+avg_sys+Total_cholesterol+HDL_cholesterol+pack_years+"
                   f"rheum_arthritis_Y+diabetes_Y+years+rank+on_pill+treated")
        cph.fit(df, duration_col='tte', event_col='event', formula=formula)
        row = summary_and_test(cph, protein, df)
        row["formula"] = formula
        row["event"] = event
        row["covar"] = "agesex"
        row["feature"] = protein
        full.append(row)

    results = pd.DataFrame(full)
    results = pd.merge(annots, results, left_on="id", right_on="feature", how="inner")
    results.drop('feature', axis=1, inplace=True)
    results.rename(columns={"Name": "name"}, inplace=True)
    results.to_csv(path + f"/medications_full_{event}_{run}.csv")

