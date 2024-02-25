# %%
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter

# %%
analysis_type = "plotly_express"
path = "results/incremental_parallel_deaths/hosp/agesex/"
df = pd.read_csv(path + "merged_results.csv")
events = df["event"].unique()
formula_basic = "age+sex+protein"
formula_full = ("age+sex+protein+avg_sys+Total_cholesterol+HDL_cholesterol+"
                "pack_years+rheum_arthritis_Y+diabetes_Y+years+rank+on_pill")
bonf_correction = 0.05/439
threshold = bonf_correction
# df = df[df["id"] == "P01019"]
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
                                        "yes" if x["p_basic"] < threshold and x["p_full"] < threshold else "no", axis=1)
    plotting = pd.concat([plotting, tmp], axis="rows")

plotting["att"] = plotting.apply(lambda row: (np.log10(row["hr_full"])*100)/np.log10(row["hr_basic"]), axis=1)
plotting.to_csv(path + "plotting_df.csv", index=False)
# %%
if analysis_type == "plotly_express":
    import plotly.express as px
    import plotly.io as pio

    pio.renderers.default = "browser"
    fig = px.scatter(plotting, x="hr_basic", y="hr_full", opacity=0.7, color="both_significant",
                     error_x="err_basic", error_y="err_full", hover_name="name_x",
                     facet_col="event", facet_col_wrap=2,
                     labels={
                         "hr_basic": "HR basic (95% CI)",
                         "hr_full": "HR full (95% CI)",
                         "event": "Event",
                         "both_significant": "Both HRs significant (Bonferroni)"
                     })
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    fig.show()
elif analysis_type == "matplotlib":
    plt.errorbar(plotting["hr_basic"], plotting["hr_full"],
                 xerr=plotting["err_basic"],
                 yerr=plotting["err_full"], fmt='o', alpha=0.6)
    plt.xlabel("HR basic")
    plt.ylabel("HR full")
    plt.title("")
    plt.show()


# Protein * sex interaction term effect
# Protein: P00751.H7C5H1.E7ETN3.B4E1Z4
#%%
df = pd.read_csv("results/cox/hosp/prepped/cox_hosp_hf_prepped.csv")
proteins = pd.read_csv('results/cox/hosp/prepped/proteins_hosp_all_events_scaled_8660.csv')
protein = "P00751.H7C5H1.E7ETN3.B4E1Z4"
df = pd.merge(df, proteins, on="id")
print(df.shape)

cph = CoxPHFitter()
formula = (f"age+sex*{protein}+avg_sys+Total_cholesterol+HDL_cholesterol+pack_years+"
           f"rheum_arthritis_Y+diabetes_Y+years+rank+on_pill")
cph.fit(df, duration_col='tte', event_col='event', formula=formula)
cph.print_summary()

plot = cph.plot_partial_effects_on_outcome(covariates=["sex", "P00751.H7C5H1.E7ETN3.B4E1Z4"],
                                           values=[["M", -3], ["M", 3],
                                                   ["F", -3], ["F", 3]], cmap='coolwarm', plot_baseline=False)
plt.ylabel("HF-free survival")
plt.xlabel("Follow up (years)")
plt.legend(["Males, Proteins 2,3 (low)", "Males, Proteins 2,3 (high)",
            "Females, Proteins 2,3 (low)", "Females, Proteins 2,3 (high)"])
# plt.savefig("plots/survival.png", dpi=600)
plt.show()
# cph.hazard_ratios_.to_csv("results/incremental_parallel/hosp/agesex_interaction/hazard_ratios.csv")
# cph.summary.to_csv("results/incremental_parallel/hosp/agesex_interaction/summary.csv")


#%% check the plotting function

# Separate data by sex
male_data_pg_low = df[(df['sex'] == 'M') & (df['P00751.H7C5H1.E7ETN3.B4E1Z4'] < -1.5)]
male_data_pg_high = df[(df['sex'] == 'M') & (df['P00751.H7C5H1.E7ETN3.B4E1Z4'] > 1.5)]
female_data_pg_low = df[(df['sex'] == 'F') & (df['P00751.H7C5H1.E7ETN3.B4E1Z4'] < -1.5)]
female_data_pg_high = df[(df['sex'] == 'F') & (df['P00751.H7C5H1.E7ETN3.B4E1Z4'] > 1.5)]

# Fit Kaplan-Meier estimator for each sex
kmf_male_low = KaplanMeierFitter()
kmf_male_high = KaplanMeierFitter()
kmf_female_low = KaplanMeierFitter()
kmf_female_high = KaplanMeierFitter()

kmf_male_low.fit(male_data_pg_low['tte'], event_observed=male_data_pg_low['event'], label='Male, PG (low)')
kmf_male_high.fit(male_data_pg_high['tte'], event_observed=male_data_pg_high['event'], label='Male, PG (high)')
kmf_female_low.fit(female_data_pg_low['tte'], event_observed=female_data_pg_low['event'], label='Female, PG (low)')
kmf_female_high.fit(female_data_pg_high['tte'], event_observed=female_data_pg_high['event'], label='Female, PG (high)')


# Plot Kaplan-Meier curves for each sex
plt.figure(figsize=(10, 6))

kmf_male_low.plot()
kmf_male_high.plot()
kmf_female_low.plot()
kmf_female_high.plot()

plt.title('Kaplan-Meier Survival Curve by Sex')
plt.xlabel('Time')
plt.ylabel('Survival Probability')
plt.grid(False)
plt.legend()
plt.show()

#%%
df = pd.read_csv("results/cox/hosp/prepped/cox_hosp_hf_prepped.csv")
proteins = pd.read_csv('results/cox/hosp/prepped/proteins_hosp_all_events_scaled_8660.csv')
protein = "P00751.H7C5H1.E7ETN3.B4E1Z4"
df = pd.merge(df, proteins, on="id")

