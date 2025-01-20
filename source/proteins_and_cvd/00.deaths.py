import pyreadr
import pandas as pd
from lib.pandas_ext import two_dfs_merge

# %%
# Read in pheno and proteins files, to get records with data in both (inner join)
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
all_CDC = ["I00", "I01", "I02", "I03", "I04", "I05", "I06", "I07", "I08", "I09", "I10",
           "I11", "I12", "I13", "I14", "I15", "I16", "I17", "I18", "I19", "I20", "I21",
           "I22", "I23", "I24", "I25", "I26", "I27", "I28", "I29", "I30", "I31", "I32",
           "I33", "I34", "I35", "I36", "I37", "I38", "I39", "I40", "I41", "I42", "I43",
           "I44", "I45", "I46", "I47", "I48", "I49", "I50", "I51", "I52", "I53", "I54",
           "I55", "I56", "I57", "I58", "I59", "I60", "I61", "I62", "I63", "I64", "I65",
           "I66", "I67", "I68", "I69", "I70", "I71", "I72", "I73", "I74", "I75", "I76",
           "I77", "I78"]

all_welsh = ["I00", "I01", "I02", "I03", "I04", "I05", "I06", "I07", "I08", "I09", "I10",
           "I11", "I12", "I13", "I14", "I15", "I16", "I17", "I18", "I19", "I20", "I21",
           "I22", "I23", "I24", "I25", "I26", "I27", "I28", "I29", "I30", "I31", "I32",
           "I33", "I34", "I35", "I36", "I37", "I38", "I39", "I40", "I41", "I42", "I43",
           "I44", "I45", "I46", "I47", "I48", "I49", "I50", "I51", "I52", "I53", "I54",
           "I55", "I56", "I57", "I58", "I59", "I60", "I61", "I62", "I63", "I64", "I65",
           "I66", "I67", "I68", "I69", "I70", "I71", "I72", "I73", "I74", "I75", "I76",
           "I77", "I78", "I79", "I80", "I81", "I90", "I91", "I92", "I93", "I94", "I95",
           "I96", "I97", "I98", "I99"]

# %%
cvd_deaths = deaths[deaths['value'].isin(heart + hypertension + cerebrovascular)]
cvd_deaths = cvd_deaths.drop_duplicates(subset='id')  # 561

# %%
cvd_deaths2 = deaths[deaths['value'].isin(all_CDC)]
cvd_deaths2 = cvd_deaths2.drop_duplicates(subset='id')  # 583

# %%
cvd_deaths3 = deaths[deaths['value'].isin(all_welsh)]
cvd_deaths3 = cvd_deaths3.drop_duplicates(subset='id')  # 583

# %%
cvd_deaths.to_csv("data/phenotypes/2025-01-13_cvd_deaths_welsh.csv", index=False)
cvd_deaths.to_csv("data/phenotypes/2024-02-23_cvd_deaths.csv", index=False)

# %%
# prepare cvd death phenotype
cvd_deaths_df = pd.merge(pheno, cvd_deaths3, on="id")  # 451 deaths in merged df vs 469 with Welsh
cvd_deaths_df = pd.DataFrame({"id": cvd_deaths_df["id"],
                              "dt1_ym": cvd_deaths_df["dod_ym_y"],
                              "gs_appt": cvd_deaths_df["gs_appt"],
                              "incident": [1] * len(cvd_deaths_df),
                              "Disease": ["CVD_death"] * len(cvd_deaths_df),
                              "Source": ["Secondary_Care"] * len(cvd_deaths_df),
                              "GP_Consent": [1] * len(cvd_deaths_df)
                              })
cvd_deaths_df = cvd_deaths_df.drop_duplicates(subset='id')  # 469 deaths in merged df
cvd_deaths_df.to_csv("data/phenotypes/2025-01-13_cvd_deaths_df_welsh.csv", index=False)

# %%
# prepare death as phenotype
# I don't have to change deaths as such, it's only the CVD death definition that's changed
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

deaths_df.to_csv("data/phenotypes/2025-01-13_deaths_df.csv", index=False)

