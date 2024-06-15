# %%
import os
import sys
import multiprocessing as mp
from tqdm import tqdm
from lifelines import CoxPHFitter

sys.path.append('/Cluster_Filespace/Marioni_Group/Ola/Code/general/projects/proteins')
from lib.cox import extract_cox_coefs, summary_and_test
import pandas as pd


def get_formulae(run="agesex", additional_params=None):
    method = '+' if run == "agesex" else '*'
    covars = ['age+sex', 'avg_sys', 'Total_cholesterol',
              'HDL_cholesterol', 'pack_years', 'rheum_arthritis_Y',
              'diabetes_Y', 'years', 'rank', 'on_pill']
    formulae = []

    for i in range(0, len(covars)):
        if i == 0:
            formulae.append(covars[i] + f"{method}protein")
        else:
            formulae.append(formulae[i - 1] + "+" + covars[i])
    return formulae


def event_dict(flag, events):
    ret_dict = {}

    for event in events:
        cox_path = f"results/cox/{flag}/prepped/cox_{flag}_{event}_prepped.csv"
        ret_dict[event] = pd.read_csv(cox_path)
        ret_dict[event].set_index("id", inplace=True)

    return ret_dict


def process_proteins(annots, events, feature, flag, interesting_events, protein, proteins, run):

    for event in interesting_events:
        path = f'results/incremental_parallel/{flag}/{run}/{event}'
        if not os.path.exists(path):
            os.makedirs(path)
            print(f"Path: {path} created!")

        formulae = get_formulae(run)

        full = []
        for formula in formulae:
            cox = events[event]
            feature = feature.replace('protein', protein)
            df = pd.merge(proteins, cox, how="inner", left_index=True, right_index=True)
            cph = CoxPHFitter()
            cph.fit(df, duration_col='tte', event_col='event', formula=formula.replace('protein', protein))
            row = summary_and_test(cph, feature, df)
            row["formula"] = formula
            row["event"] = event
            row["covar"] = feature
            row["feature"] = protein
            full.append(row)

        results = pd.DataFrame(full)
        results = pd.merge(annots, results, left_on="id", right_on="feature", how="inner")
        results.drop('feature', axis=1, inplace=True)
        results.rename(columns={"Name": "name"}, inplace=True)
        results.to_csv(path + f"/{protein}.csv")


def main():

    flag = "hosp"  # hosp_gp, hosp, hosp_gp_cons
    run = "agesex_interaction"
    # feature = "sex[T.M]:protein"
    #run = "agesex"
    feature = "protein"
    cores = int(mp.cpu_count() * 0.8)
    print(f"Used cores: {cores}")
    proteins = pd.read_csv('results/cox/hosp/prepped/proteins_hosp_all_events_scaled_8491.csv')
    annots = pd.read_csv("data/annotations/short_annots.csv")
    proteins.set_index("id", inplace=True)
    interesting_events = ["myocardial_infarction", "isch_stroke", "hf", "chd_nos",
                           "tia", "composite_CVD", "CVD_death", "death"]
    events = event_dict(flag, interesting_events)

    path = f'results/incremental_parallel/{flag}/{run}'
    if not os.path.exists(path):
        os.makedirs(path)
        print(f"Path: {path} created!")

    pbar = tqdm(total=len(proteins.columns))

    def update(*a):
        pbar.update()

    with mp.Pool(cores) as pool:
        for protein in proteins.columns:
            pool.apply_async(process_proteins,
                             args=(annots, events, feature, flag, interesting_events, protein, proteins, run,),
                             callback=update)
        pool.close()
        pool.join()


if __name__ == "__main__":
    main()
