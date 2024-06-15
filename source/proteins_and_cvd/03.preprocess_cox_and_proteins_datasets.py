#%%
import os
import pyreadr
import pandas as pd
import matplotlib.pyplot as plt
from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test
from lib.stats import scale, summary
from lib.pandas_ext import two_dfs_merge
from lib.parquet_helper import write_parquet
from lib.outliers import outlier_trim
import numpy as np


flag = "hosp"  # hosp_gp, hosp, hosp_gp_cons
proteins = pyreadr.read_r("data/phenotypes/GS_ProteinGroups_RankTransformed_23Aug2023.rds")[None]
interesting_events = ["myocardial_infarction", "isch_stroke", "hf", "chd_nos",
                      "tia", "composite_CVD", "CVD_death", "death"]

# initial subsetting
event = "hf"
cox_path = f"results/cox/{flag}/{flag}_{event}.csv"
cox = pd.read_csv(cox_path)
cox = cox[['id', 'age', 'sex']]  # for readability
cox.set_index("id", inplace=True)
proteins.set_index("id", inplace=True)

cox, proteins = two_dfs_merge(cox, proteins)  # both 13374 records, need to recalculate events
indexes = cox.index

path = f'results/cox/{flag}/prepped/'
if not os.path.exists(path):
    os.makedirs(path)
    print(f"Path: {path} created!")

#%%
# prep for cvd ~ age + sex + protein + risk_factors + on_pill
for event in interesting_events:
    cox_path = f"results/cox/{flag}/{flag}_{event}.csv"
    cox = pd.read_csv(cox_path)
    cox = cox[cox["id"].isin(indexes)]
    cox = cox[['id', 'age', 'sex', 'avg_sys', 'Total_cholesterol',
               'HDL_cholesterol', 'pack_years', 'rheum_arthritis_Y', 'diabetes_Y',
               'bmi', 'years', 'rank', 'event', 'tte', 'on_pill']]
    cox["sex"] = pd.Categorical(cox["sex"])
    cox_plotting = cox.copy()

    cox["HDL_cholesterol"] = outlier_trim(cox["HDL_cholesterol"], cut=4)
    cox["Total_cholesterol"] = outlier_trim(cox["Total_cholesterol"], cut=4)
    cox['bmi'] = np.where((cox['bmi'] < 18) | (cox['bmi'] > 50), np.NaN, cox['bmi'])
    cox['avg_sys'] = outlier_trim(cox["avg_sys"], cut=4)
    cox['pack_years'] = outlier_trim(cox["pack_years"], cut=4)
    cox.dropna(inplace=True)

    cox = cox.query('age <= 69 and age >= 39')
    cph = CoxPHFitter()
    cph.fit(cox, duration_col='tte', event_col='event',
            formula=f"age+sex+Total_cholesterol+HDL_cholesterol+"
                    f"avg_sys+pack_years+rheum_arthritis_Y+diabetes_Y+"
                    f"rank+on_pill")
    test = proportional_hazard_test(cph, cox, time_transform="km")
    print(event)
    print(f"Number of events {np.sum(cox.event)}")
    print(f"Length: {len(cox)}")
    print(test.summary)
    print("====\n\n")

    cox.to_csv(f"results/cox/hosp/prepped/cox_{flag}_{event}_prepped.csv", index=False)

#%%
cox.set_index("id", inplace=True)
cox, proteins = two_dfs_merge(cox, proteins)  # both 8660 records
proteins = scale(proteins)
proteins.to_csv(f"results/cox/hosp/prepped/proteins_{flag}_all_events_scaled_8491.csv")


# basic models solved!
#%%
plt.hist(cox["bmi"])
plt.title("BMI")
plt.show()

plt.hist(cox["avg_sys"])
plt.title("AVG_SYS")
plt.show()

plt.hist(cox["pack_years"])
plt.title("PY")
plt.show()

plt.boxplot(cox_plotting["Total_cholesterol"])
plt.title("Total cholesterol")
plt.show()

plt.boxplot(cox_plotting["HDL_cholesterol"])
plt.title("HDL cholesterol")
plt.show()
#%%
cox_plotting["Total_cholesterol"] = outlier_trim(np.log10(cox_plotting["Total_cholesterol"]), cut=3)
cox_plotting["HDL_cholesterol"] = outlier_trim(np.log10(cox_plotting["HDL_cholesterol"]), cut=3)
cox_plotting.dropna(inplace=True)

plt.boxplot(cox_plotting["HDL_cholesterol"])
plt.title("HDL cholesterol")
plt.show()

plt.boxplot(cox_plotting["Total_cholesterol"])
plt.title("Total cholesterol")
plt.show()
