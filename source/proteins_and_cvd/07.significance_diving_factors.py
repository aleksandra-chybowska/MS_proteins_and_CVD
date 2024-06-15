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


# %%
# we have an error here - not age * sex, it should be protein * sex
def get_formulae():
    basic = "age+sex+protein"
    covars = ['avg_sys', 'Total_cholesterol',
              'HDL_cholesterol', 'pack_years', 'rheum_arthritis_Y',
              'diabetes_Y', 'years', 'rank', 'on_pill']
    formulae = [basic]

    for i in range(0, len(covars)):
        formulae.append(basic + f"+{covars[i]}")
    return formulae


def main():
    # %%
    flag = "hosp"  # hosp_gp, hosp, hosp_gp_cons
    run = "agesex"
    input_path = "results/incremental_parallel_correct/hosp/agesex/"

    proteins = pd.read_csv('results/cox/hosp/prepped/proteins_hosp_all_events_scaled_8660.csv')
    annots = pd.read_csv("data/annotations/short_annots.csv")
    proteins.set_index("id", inplace=True)
    plotting = pd.read_csv(input_path + "plotting_df.csv")
    interesting_events = plotting["event"].unique()

    path = f'results/attenuation_factor/{run}/{flag}'
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
        formulae = get_formulae()

        sig_proteins = plotting[plotting["event"] == event]
        sig_proteins = sig_proteins[sig_proteins["both_significant"] == "no"]
        ids = sig_proteins["id"]

        full = []

        for protein in ids:
            for formula in formulae:
                print(protein)
                cph = CoxPHFitter()
                cph.fit(df, duration_col='tte', event_col='event', formula=formula.replace('protein', protein))
                row = summary_and_test(cph, protein, df)
                row["formula"] = formula
                row["event"] = event
                row["covar"] = protein
                row["feature"] = protein
                full.append(row)
                # concordance?

        results = pd.DataFrame(full)
        results = pd.merge(annots, results, left_on="id", right_on="feature", how="inner")
        results.drop('feature', axis=1, inplace=True)
        results.rename(columns={"Name": "name"}, inplace=True)
        results.to_csv(path + f"/attenuation_factor_{event}_{run}.csv")


if __name__ == "__main__":
    main()
