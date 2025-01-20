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
proteins = pd.read_csv('results/cox/40-69/proteins_hosp_all_events_scaled_8343.csv')
# sanitized_proteins = {col: col.replace('.', '_') for col in proteins.columns}
# proteins = proteins.rename(columns=sanitized_proteins)
proteins.set_index("id", inplace=True)

sig = pd.read_csv("results/incremental_parallel/hosp/agesex/40-69/merged_results_bonf_significant_full.csv")
#sig['id'] = sig['id'].str.replace('.', '_')
proteins = proteins.loc[:, sig["id"].unique()]
annots = pd.read_csv("data/annotations/short_annots.csv")

path = f'results/medications/{run}/{flag}'
if not os.path.exists(path):
    os.makedirs(path)
    print(f"Path: {path} created!")

# look only at events for which significant associations were reported
events = sig["event"].unique()
# %%
for event in events:
    print(event)
    cox_path = f"results/cox/40-69/cox_{flag}_{event}_prepped.csv"
    cox = pd.read_csv(cox_path)
    cox["treated"] = 0
    cox.loc[cox["id"].isin(treated), "treated"] = 1
    cox.set_index("id", inplace=True)
    # 1786 treated of n = 8660, 1761 in m=8343
    df = pd.merge(proteins, cox, how="inner", left_index=True, right_index=True)
    ids = proteins.columns
    full = []

    for item in ids:
        print(item)
        protein = item
        cph = CoxPHFitter()
        formula = ("age+sex+protein+avg_sys+Total_cholesterol+HDL_cholesterol+pack_years+rheum_arthritis_Y+diabetes_Y+years+rank+on_pill+treated")
        cph.fit(df, duration_col='tte', event_col='event', formula=formula.replace('protein', protein))
        row = summary_and_test(cph, item, df)
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

