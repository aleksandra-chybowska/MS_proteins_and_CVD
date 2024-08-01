import numpy as np
import pandas as pd


# %%
path = "results/incremental_parallel/hosp/agesex/40-69/"
df = pd.read_csv(path + "merged_results.csv")
events = df["event"].unique()
formula_basic = "age+sex+protein"
formula_full = ("age+sex+protein+avg_sys+Total_cholesterol+HDL_cholesterol+"
                "pack_years+rheum_arthritis_Y+diabetes_Y+years+rank+on_pill")
bonf_correction = 0.05 / 439
threshold = bonf_correction
plotting = pd.DataFrame()
# %%
for event in events:
    tmp = df[df["event"] == event]
    basic = tmp[tmp["formula"] == formula_basic]
    basic = basic[["id", "name", "hr", "lci", "uci", "p"]]
    basic.rename(columns={"hr": "hr_basic", "lci": "lci_basic", "uci": "uci_basic", "p": "p_basic"},
                 inplace=True)
    basic["err_basic"] = (basic["uci_basic"] - basic["lci_basic"]) / 2
    significant_basic = basic[basic["p_basic"] < threshold]
    names = set(significant_basic["name"])

    full = tmp[tmp["formula"] == formula_full]
    full = full[["id", "hr", "lci", "uci", "event", "name", "p"]]
    full.rename(columns={"hr": "hr_full", "lci": "lci_full", "uci": "uci_full", "p": "p_full"},
                inplace=True)
    full["err_full"] = (full["uci_full"] - full["lci_full"]) / 2
    significant_full = full[full["p_full"] < threshold]
    names.update(set(significant_full["name"]))

    tmp = pd.merge(basic, full, on="id", how="inner")
    tmp = tmp.query("name_x in @names").copy()
    tmp["both_significant"] = tmp.apply(lambda x:
                                        "yes" if x["p_basic"] < threshold and x["p_full"] < threshold else "no",
                                        axis=1)
    plotting = pd.concat([plotting, tmp], axis="rows")

plotting["att"] = plotting.apply(lambda row: 100 - ((np.log(row["hr_full"]) * 100) / np.log(row["hr_basic"])),
                                 axis=1)
plotting.to_csv(path + "plotting_df_new.csv", index=False)
#
# for event in events:
#     sig_proteins = plotting[plotting["event"] == event]
#     sig_proteins = sig_proteins[sig_proteins["both_significant"] == "no"]
#     ids = sig_proteins["id"]
#     tmp = df[(df["event"] == event) & (df["id"].isin(ids))]
#
#     for id in ids:
#         significant = "yes"
#         # get all covars as table
#         # run cox model, basic
#         # save results
#         # run cox models with first covar added
#         # check if protein is significant
#         # if significant, continue
#         # if not, save the protein to a dataframe - show covar. Basic model + results with that group added.
