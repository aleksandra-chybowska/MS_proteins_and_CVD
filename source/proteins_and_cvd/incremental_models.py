# This script compares the following coxPH models:
# a)	Null models, adjusted for age, sex and an individual protein abundance.
# b)	Full models, that will include CVD risk factors (such as age, sex, systolic blood pressure,
#       total cholesterol, HDL-cholesterol, LDL-cholesterol, pack years of smoking, family history of CVD,
#       rheumatoid arthritis, diabetes, BMI, years of education, SIMD score) and individual protein abundance.

# %%
import os
import sys
import pandas as pd
from lifelines import CoxPHFitter

sys.path.append('/Cluster_Filespace/Marioni_Group/Ola/Code/general/projects/proteins')
from lib.cox import extract_cox_coefs, summary_and_test


def get_formulae(run="agesex", additional_params=None):

    first_covar = "age+sex"
    if run == "agesex_interaction":
        first_covar = "age*sex"

    covars = [first_covar, 'avg_sys', 'Total_cholesterol',
              'HDL_cholesterol', 'pack_years', 'rheum_arthritis_Y',
              'diabetes_Y', 'years', 'rank', 'on_pill']
    formulae = []

    for i in range(0, len(covars)):
        if i == 0:
            formulae.append(covars[i])
        else:
            formulae.append(formulae[i - 1] + "+" + covars[i])

    return formulae


def main():
    # %%
    flag = "hosp"  # hosp_gp, hosp, hosp_gp_cons
    run = "agesex_interaction"

    proteins = pd.read_csv('results/cox/hosp/prepped/proteins_hosp_all_events_scaled_8660.csv')
    annots = pd.read_csv("data/annotations/short_annots.csv")
    proteins.set_index("id", inplace=True)
    interesting_events = ["myocardial_infarction", "isch_stroke", "hf", "chd_nos",
                          "tia", "composite_CVD", "CVD_death"]

    path = f'results/incremental_models/{run}/{flag}'
    if not os.path.exists(path):
        os.makedirs(path)
        print(f"Path: {path} created!")

    # %%
    for event in interesting_events:
        print(event)
        print()
        cox_path = f"results/cox/{flag}/prepped/cox_{flag}_{event}_prepped.csv"
        cox = pd.read_csv(cox_path)
        cox.set_index("id", inplace=True)
        df = pd.merge(proteins, cox, how="inner", left_index=True, right_index=True)

        # Here implement all these fantastic types of models :)
        # age + sex + Total_chol...
        # age + sex + Total_chol... + HDL_chol...
        formulae = get_formulae(run)

        for formula in formulae:
            path = f'results/incremental_models/{run}/{flag}'
            tmp = formula.replace('*', 'x')
            path = path + "/" + tmp.replace('+', '_') + "/"

            if not os.path.exists(path):
                os.makedirs(path)
                print(f"Path: {path} created!")

            full = []
            for protein in proteins.columns:
                print(protein)
                cph = CoxPHFitter()
                cph.fit(df, duration_col='tte', event_col='event', formula=formula + f"+{protein}")
                row = summary_and_test(cph, protein, df)
                full.append(row)
                # concordance?
        # %%
            results = pd.DataFrame(full)
            results = pd.merge(annots, results, left_on="id", right_on="feature", how="inner")
            results.drop('feature', axis=1, inplace=True)
            results.rename(columns={"Name": "name"}, inplace=True)
            results.to_csv(path + f"/full_{event}_{run}.csv")


if __name__ == "__main__":
    main()
