#%%
import pyreadr
import pandas as pd
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test
from lib.stats import scale, summary
from lib.pandas_ext import two_dfs_merge
from lib.parquet_helper import write_parquet, read_parquet
import numpy as np
from sklearn.preprocessing import StandardScaler

flag = "hosp"  # hosp_gp, hosp, hosp_gp_cons
proteins = pyreadr.read_r("data/phenotypes/GS_ProteinGroups_RankTransformed_23Aug2023.rds")[None]
interesting_events = ["myocardial_infarction", "isch_stroke", "hf", "chd_nos",
                      "tia", "composite_CVD", "CVD_death"]

# initial subsetting
event = "hf"
cox_path = f"results/cox/{flag}/{flag}_{event}.csv"
cox = pd.read_csv(cox_path)
cox = cox[['id', 'age', 'sex']]  # for readability
cox.set_index("id", inplace=True)
proteins.set_index("id", inplace=True)

cox, proteins = two_dfs_merge(cox, proteins)  # both 13374 records, need to recalculate events
indexes = cox.index
#%%
# first, lets see if the basic models meet cox assumptions and if not, lets fix them
for event in interesting_events:
    cox_path = f"results/cox/{flag}/{flag}_{event}.csv"
    cox = pd.read_csv(cox_path)
    cox = cox[cox["id"].isin(indexes)]
    cox = cox[['id', 'age', 'sex', 'avg_sys', 'Total_cholesterol',
               'HDL_cholesterol', 'pack_years', 'rheum_arthritis_Y', 'diabetes_Y',
               'bmi', 'years', 'rank', 'event', 'tte', 'on_pill']]
    cox["sex"] = pd.Categorical(cox["sex"])
    cox = cox[cox["HDL_cholesterol"] < 3.5]  # put here 3.0 if you want a stable model
    cox = cox[cox["Total_cholesterol"] < 8.0]  # put here 6.0 if you want a stable model

    cox["HDL_cholesterol"] = np.log10(cox["HDL_cholesterol"])
    cox["Total_cholesterol"] = np.log10(cox["Total_cholesterol"])

    # lets start from a basic cox model
    # cox = cox.query('age <= 60 and age >= 40')
    cox = cox.query('age <= 69 and age >= 39')
    cph = CoxPHFitter()
    cph.fit(cox, duration_col='tte', event_col='event',
            formula=f"age+sex+Total_cholesterol+HDL_cholesterol+"
                    f"avg_sys+pack_years+rheum_arthritis_Y+diabetes_Y+"
                    f"rank+on_pill")
    test = proportional_hazard_test(cph, cox, time_transform="km")
    print(event)
    print(f"Number of events {np.sum(cox.event)}")
    print(test.summary)
    print("====\n\n")


# basic models solved!
#%%
# plt.hist(cox["bmi"])
# plt.show()

event = "hf"  # any event would do, these ds's contain the same individuals


#### correct protein dataset here based on filtering in the loop

# serious trust issues are evident
common = pd.merge(cox, proteins, left_index=True, right_index=True)
np.allclose(proteins.iloc[:, 0], common.iloc[:, 2])
(proteins.iloc[:, 0] == common.iloc[:, 2]).all()

scaled_proteins = scale(proteins)
# very evident indeed
summary(scaled_proteins.iloc[:, 0:4])

scaled_proteins.reset_index(inplace=True) # I think it needs "inplace". I left the index as it was for now
write_parquet(scaled_proteins, "data/transformed_input/cox_analysis_proteins_scaled_13374.parquet")
