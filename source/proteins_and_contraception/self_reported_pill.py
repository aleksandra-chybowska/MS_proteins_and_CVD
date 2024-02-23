# %%
import pyreadr
import pandas as pd
import numpy as np
import statsmodels.api as sm
from lib.lm import results_summary_to_dataframe
from helpers.womans_health import get_on_pill
from lib.parquet_helper import read_parquet


pheno = pyreadr.read_r("data/phenotypes/GS_phenos_internal_with_DST_28Nov2023_REM.rds")[None]
proteins = pyreadr.read_r("data/phenotypes/GS_ProteinGroups_RankTransformed_23Aug2023.rds")[None]
annots = pd.read_csv("data/annotations/short_annots.csv")
pill = pd.read_table("data/disease/womans_phenotypes/GS_womens_phenotypes_v2v5combined.txt")

# %%
cols = ["id", "age", "sex", "taken_cont", "age_started_cont", "years_taking_cont"]
females = read_parquet("data/transformed_input/females_and_proteins_scaled.parquet")

# %%
pd.crosstab(columns="count", index=females["on_pill"])  # 2082 females that could be on pill

cols.append("on_pill")
females_pheno = females[cols]
females_proteins = females[females.columns.difference(cols)]

#%%
df = pd.DataFrame()
for protein in females_proteins.columns:

    # define response variable
    y = females_proteins[[protein]]
    x = sm.add_constant(females_pheno[["on_pill", "age"]])
    feature = "on_pill"
    # fit linear regression model
    model = sm.OLS(y, x).fit()
    summary = results_summary_to_dataframe(model)
    new_data = pd.DataFrame({"Protein": [protein],
                             "Coef": [summary.loc[feature, 'coeff']],
                             "P": [summary.loc[feature, 'pvals']]})
    df = pd.concat([df, new_data], ignore_index=True)

df = pd.merge(df, annots, left_on="Protein", right_on="id", right_index=False).drop(columns='id')
# %%


df = df.sort_values(by='P', ascending=True)
df.loc[df.P < 0.05]  # 204 rows

bonf = 0.05/len(females_proteins)  # 66 rows
corrected = df.loc[df.P < bonf]

df.to_csv("results/proteins_that_vary_with_pill.csv", index=False)
corrected.to_csv("results/proteins_that_vary_with_pill_bonf_corrected.csv", index=False)
