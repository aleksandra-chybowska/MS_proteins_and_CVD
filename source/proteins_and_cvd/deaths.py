import pyreadr
import pandas as pd
from lib.pandas_ext import two_dfs_merge

pheno = pd.read_csv("data/transformed_input/generic_pheno.csv")  # based on deaths_old, there's 1122 deaths in this file
proteins = pyreadr.read_r("data/phenotypes/GS_ProteinGroups_RankTransformed_23Aug2023.rds")[None]
pheno.set_index("id", inplace=True)
proteins.set_index("id", inplace=True)

pheno, proteins = two_dfs_merge(pheno, proteins)  # both 13374 records, need to recalculate events
pheno = pheno.reset_index()
proteins = proteins.reset_index()
# %%
deaths = pyreadr.read_r("data/phenotypes/2024-02-19_proteomics_mortality_phenodata_for_ola.rds")[None]
deaths.set_index("id", inplace=True)
pivot = deaths.loc[:, "cause":"sec9"]
pivot = pivot.reset_index()
pivot = pd.melt(pivot, id_vars=["id"], var_name="column", value_name="value")
pivot.drop(columns=["column"], inplace=True)
pivot.dropna(inplace=True)

# %%
deaths = deaths.reset_index()
deaths = deaths[['id', 'dod_ym', 'dead']]
deaths.dropna(inplace=True)  # 1245 new deaths
pheno['dead'].value_counts()  # 870 old deaths
deaths = pd.merge(deaths, pivot, on="id", how="left")
deaths['value'] = deaths['value'].str.extract(r'(...)')  # substring entire column, 3 letters

deaths.to_csv("data/phenotypes/2024-02-23_deaths.csv", index=False)

# %%
# CVD deaths
heart = ["I00", "I01", "I02", "I03", "I04", "I05", "I06", "I07", "I08", "I09", "I11",
         "I13", "I20", "I21", "I22", "I23", "I24", "I25", "I26", "I27", "I28", "I29",
         "I30", "I31", "I32", "I33", "I34", "I35", "I36", "I37", "I38", "I39", "I40",
         "I41", "I42", "I43", "I44", "I45", "I46", "I47", "I48", "I49", "I50", "I51"]
hypertension = ["I10", "I12", "I15"]
cerebrovascular = ["I60", "I61", "I62", "I63", "I64", "I65", "I66", "I67", "I68", "I69"]
# %%
cvd_deaths = deaths[deaths['value'].isin(heart + hypertension + cerebrovascular)]
cvd_deaths = cvd_deaths.drop_duplicates(subset='id')  # 561
cvd_deaths.to_csv("data/phenotypes/2024-02-23_cvd_deaths.csv", index=False)

# check with previous file
# %%
pepsi = pd.read_csv("data/disease/cvd_deaths.csv")
pepsi = pepsi[pepsi['value'].isin(heart + hypertension + cerebrovascular)]
pepsi = pepsi.drop_duplicates(subset='id')  # 693
pepsi = pepsi[pepsi['id'].isin(pheno['id'])]
pepsi = pepsi.drop(columns=['value'])  # 365 with proteins measured
check = pd.merge(pepsi, cvd_deaths, how='outer', on='id')  # NAs are from pepsi - all is fine

# %%
# prepare cvd death phenotype
cvd_deaths_df = pd.merge(pheno, cvd_deaths, on="id")  # 451 deaths in merged df
cvd_deaths_df = pd.DataFrame({"id": cvd_deaths_df["id"],
                              "dt1_ym": cvd_deaths_df["dod_ym_y"],
                              "gs_appt": cvd_deaths_df["gs_appt"],
                              "incident": [1] * len(cvd_deaths_df),
                              "Disease": ["CVD_death"] * len(cvd_deaths_df),
                              "Source": ["Secondary_Care"] * len(cvd_deaths_df),
                              "GP_Consent": [1] * len(cvd_deaths_df)
                              })
cvd_deaths_df = cvd_deaths_df.drop_duplicates(subset='id')  # 451 deaths in merged df
cvd_deaths_df.to_csv("data/phenotypes/2024-02-23_cvd_deaths_df.csv", index=False)

# %%
# prepare death as phenotype
deaths_df = deaths.drop_duplicates(subset='id')  # 1245
deaths_df = pd.merge(pheno, deaths_df, on="id")  # 1022
deaths_df = pd.DataFrame({"id": deaths_df["id"],
                          "dt1_ym": deaths_df["dod_ym_y"],
                          "gs_appt": deaths_df["gs_appt"],
                          "incident": [1] * len(deaths_df),
                          "Disease": ["death"] * len(deaths_df),
                          "Source": ["Secondary_Care"] * len(deaths_df),
                          "GP_Consent": [1] * len(deaths_df)
                          })

deaths_df.to_csv("data/phenotypes/2024-02-23_deaths_df.csv", index=False)

