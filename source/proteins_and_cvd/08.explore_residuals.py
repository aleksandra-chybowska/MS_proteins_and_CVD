import os

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from lib.cox import extract_cox_coefs, summary_and_test

# Set working directory
os.chdir("C:/Users/s1654019/Projects/python/proteins")

flag = "hosp"
run = "agesex"
type = "40-69"
bonf = 0.05 / 439

# %%
path = f"results/incremental_parallel/{flag}/{run}/{type}/"
ds = pd.read_csv(f"{path}merged_results.csv")

full_mod = ("age+sex+protein+avg_sys+Total_cholesterol+HDL_cholesterol+"
            "pack_years+rheum_arthritis_Y+diabetes_Y+years+rank+on_pill")
sig_proteins = ds.query(f'formula == "{full_mod}" and p < {bonf}')

proteins = pd.read_csv('results/cox/40-69/proteins_hosp_all_events_scaled_8332.csv')
cols = ['id'] + sig_proteins['id'].tolist()
proteins = proteins[proteins.columns.intersection(cols)]
results = pd.DataFrame()

# for i, row in sig_proteins.iterrows():
row = sig_proteins.iloc[0]
protein_name = row['id']
protein_df = proteins[['id', protein_name]].rename(columns={protein_name: 'protein'})

event = row['event']
cox_path = f"results/cox/{type}/cox_{flag}_{event}_prepped.csv"
cox = pd.read_csv(cox_path)
cox['sex'] = cox['sex'].replace({'M': 1, 'F': 0})
cox = pd.merge(cox, protein_df, on="id")

model = CoxPHFitter()
model.fit(cox, duration_col='tte', event_col='event', formula=full_mod)

summary = model.summary
martingale_residuals = model.compute_residuals(cox, 'martingale')
r = model.compute_residuals(cox, 'deviance')
plot = r.plot.scatter(
    x='tte', y='deviance', c=np.where(r['event'], '#008fd5', '#fc4f30'),
    alpha=0.75
)
plot.figure.show()
import matplotlib.pyplot as plt

# Scatter plot of martingale residuals vs duration
plt.scatter(martingale_residuals['tte'], martingale_residuals['martingale'])
plt.xlabel('Duration')
plt.ylabel('Martingale Residuals')
plt.title('Martingale Residuals vs Duration')
plt.axhline(y=0, color='r', linestyle='--')
plt.show()
