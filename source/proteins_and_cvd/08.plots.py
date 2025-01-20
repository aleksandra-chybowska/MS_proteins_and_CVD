# %%
import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lib.cox import plot_partial_effects
import matplotlib.pyplot as plt
from contextlib import redirect_stdout


# %%
# Protein * sex interaction term effect
# Protein: P00751.H7C5H1.E7ETN3.B4E1Z4
# %%
df_hf = pd.read_csv("results/cox/40-69/cox_hosp_hf_prepped.csv")
df_death = pd.read_csv("results/cox/40-69/cox_hosp_death_prepped.csv")
proteins = pd.read_csv('results/cox/40-69/proteins_hosp_all_events_scaled_8343.csv')
annots = pd.read_csv("data/annotations/short_annots.csv")
sex_effects = pd.read_csv("results/incremental_parallel/hosp/" +
                          "agesex_interaction/40-69/merged_results_bonf_significant_full.csv")
proteins_sex_effect = sex_effects['id']

# %%
hf = pd.merge(df_hf, proteins, on="id")
death = pd.merge(df_death, proteins, on="id")

for protein in proteins_sex_effect:
    df = hf
    cph = CoxPHFitter()
    formula = (f"age+sex*{protein}+avg_sys+Total_cholesterol+HDL_cholesterol+pack_years+"
               f"rheum_arthritis_Y+diabetes_Y+years+rank+on_pill")

    cph.fit(df, duration_col='tte', event_col='event', formula=formula)
    cph.print_summary()

    covars = ["sex", protein]
    typical_range = [["M", -3], ["M", 3], ["F", -3], ["F", 3]]
    annot = annots.loc[annots["id"] == protein, "Name"].values[0]

    legend = [f"Males, {annot} (low)",
              f"Males, {annot} (high)",
              f"Females, {annot} (low)",
              f"Females, {annot} (high)"]

    plot_partial_effects(cph, covars, typical_range, legend, save=True)

    # now divide into males and females and re-run

    df_males = df[df.sex == "M"]

    cph = CoxPHFitter()
    formula = (f"age+sex+{protein}+avg_sys+Total_cholesterol+HDL_cholesterol+pack_years+"
               f"rheum_arthritis_Y+diabetes_Y+years+rank")

    cph.fit(df_males, duration_col='tte', event_col='event', formula=formula)

    with open("results/incremental_parallel/hosp/agesex_interaction/40-69/modelsummary_males.txt", 'w') as f:
        with redirect_stdout(f):
            cph.print_summary()

    df_females = df[df.sex == "F"]

    cph = CoxPHFitter()
    formula = (f"age+sex+{protein}+avg_sys+Total_cholesterol+HDL_cholesterol+pack_years+"
               f"rheum_arthritis_Y+diabetes_Y+years+rank+on_pill")

    cph.fit(df_females, duration_col='tte', event_col='event', formula=formula)

    with open("results/incremental_parallel/hosp/agesex_interaction/40-69/modelsummary_females.txt", 'w') as f:
        with redirect_stdout(f):
            cph.print_summary()



# cph.hazard_ratios_.to_csv("results/incremental_parallel/hosp/agesex_interaction/hazard_ratios.csv")
# cph.summary.to_csv("results/incremental_parallel/hosp/agesex_interaction/summary.csv")
# %% check the plotting function

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

# %%

# cross plot
analysis_type = "plotly_express"
path = "results/incremental_parallel/hosp/agesex/40-69/"
plotting = pd.read_csv(path + "plotting_df_new.csv")

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
