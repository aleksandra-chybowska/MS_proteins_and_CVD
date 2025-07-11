# %%
import os
import pyreadr
import pandas as pd
import lib.string_date as sd
from lib.cox import get_time_to_event
from lib.parquet_helper import read_parquet
# %%
proteins = pyreadr.read_r("data/phenotypes/GS_ProteinGroups_RankTransformed_23Aug2023.rds")[None]
pill = pd.read_table("data/disease/womans_phenotypes/GS_womens_phenotypes_v2v5combined.txt")
dictionary = pd.read_table("data/disease/womans_phenotypes/GS_womens_phenotypes_datadictionary.txt")
medications = pd.read_table("data/disease/medications/all_CVD_prescriptions.txt")
pheno = pyreadr.read_r("data/phenotypes/GS_phenos_internal_with_DST_28Nov2023_REM.rds")[None]  # 24079
diseases = pd.read_csv("data/phenotypes/2023-08-22_disease_codes_combined.csv")
females = read_parquet("data/transformed_input/females_on_pill.parquet")
deaths = pd.read_csv("data/phenotypes/2025-01-13_deaths_df.csv")
cvd_deaths = pd.read_csv("data/phenotypes/2025-01-13_cvd_deaths_df_welsh.csv")

# %%
flag = "hosp"  # hosp, hosp_gp_cons, hosp_gp
output = f"results/cox/basic/"

# ### phenos ###
# correlation between covariates checked in covar_corr_and_deaths.py - I will additionally check VIF too
pheno["gs_appt"] = pheno.apply(lambda row: sd.year_month_to_date(row["y_appt"], row["m_appt"]), axis="columns")
pheno = pheno[["id", "gs_appt", "age", "sex", "avg_sys", "Total_cholesterol",
               "HDL_cholesterol", "pack_years", "rheum_arthritis_Y", "diabetes_Y",
               "bmi", "years", "rank"]]
pheno = pheno.dropna()  # 17529

# %%
females = females[["id", "on_pill"]]  # 13419
pheno = pd.merge(pheno, females, on="id", how="left")
pheno["on_pill"].value_counts(dropna=False)  # 7669 NaNs
pheno.loc[pheno["sex"] == "M", "on_pill"] = 0
pheno["on_pill"].value_counts(dropna=False)  # 379 NaNs - females that had missing info for contraception
pd.crosstab(index=pheno["sex"], columns=pheno["on_pill"])
pheno = pheno.dropna()  # 17150

# %%
diseases["Disease"] = diseases["Disease"].str.lower()

score2 = pd.read_csv("data/scores/GS_score2_export.csv")
score2.value_counts(score2['notes'], dropna=False)
len(score2)  # 16934
len(score2[~score2["notes"].isna()])
score2 = score2[~score2["notes"].isna()]
score2 = score2.iloc[:, :-1]  # remove last column
score2.shape  # dim

assign = pd.read_csv("data/scores/GS_assign_export.csv")
len(assign)  # 16366
assign.value_counts(assign['notes'], dropna=False)
assign = assign[~assign["notes"].isna()]
assign = assign.iloc[:, :-1]  # remove last column
assign.shape  # 3088, 2
# %%
pheno = pd.merge(pheno, assign, on="id", how="left")
pheno = pd.merge(pheno, score2, on="id", how="left")
pheno = pd.merge(pheno, deaths[["id", "dt1_ym"]], on="id", how="left")
pheno = pheno.rename(columns={"dt1_ym": "dod_ym"})
pheno.describe()
pheno.info()
pheno["dead"] = 0
pheno.loc[~pheno['dod_ym'].isna(), 'dead'] = 1  # 1022 deaths
pheno.to_csv("data/transformed_input/generic_pheno.csv", index=False)
# %%
# diseases #
interesting_events = ["myocardial_infarction", "isch_stroke", "hf", "chd_nos", "tia"]

# cvd.query('not (Disease == "hf" and Source == "Primary_Care")')
cvd = diseases[diseases["Disease"].isin(interesting_events)]  # 5917 events
cvd = cvd.loc[~((cvd["Disease"] == "hf") & (cvd["Source"] == "Primary_Care"))]  # 5539
cvd.loc[cvd["Disease"] == "hf", "Source"].value_counts()
cvd = pd.concat([cvd, cvd_deaths], ignore_index=True) # 6008
cvd = pd.concat([cvd, deaths], ignore_index=True) # 7030
composite = cvd.query('Disease in ["isch_stroke", "chd_nos", "myocardial_infarction", "CVD_death"]').copy()
composite["Disease"] = "composite_CVD"  # 4776
# %%

cvd = pd.concat([cvd, composite], ignore_index=True)
cvd = cvd.query('incident == 1')  # 6933 + 1022 deaths = 7955

if flag == "hosp":
    cvd = cvd.query('Source == "Secondary_Care"') # 5634
if flag == "hosp_gp_cons":
    cvd = cvd.query('GP_Consent == 1')
# %%
pd.crosstab(cvd.Source, cvd.Disease, normalize="columns")
cvd['Disease'].value_counts()

dis_list = cvd.Disease.unique()
# %%
if not os.path.exists(output):
    os.makedirs(output)
    print(f"Path: {output} created!")

# %%
for outcome in dis_list:
    dis1 = cvd.query('Disease == @outcome').copy().reset_index(drop=True)  # for composite_CVD it is 2221
    dis1 = dis1.sort_values(by=["id", "dt1_ym"])
    dis1 = dis1.drop_duplicates(subset='id') # for composite CVD it is 1599
    print(f"Disease: {outcome}; N={len(dis1)}")

    out = pd.merge(pheno, dis1[["id", "dt1_ym", "Disease"]], on="id", how="left")
    out["event"] = (~out["dt1_ym"].isna()).astype(int)
    out["tte"] = out.apply(lambda row:
                           get_time_to_event(date_baseline=row["gs_appt"],
                                             date_event=row["dt1_ym"],
                                             date_censor="202308",
                                             date_death=row["dod_ym"]), axis="columns")
    out["age_at_event"] = out.age + out.tte
    out.to_csv(f"{output}/{flag}_{outcome}.csv", na_rep='NA', index=False)

