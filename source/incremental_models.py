# This script compares:
# a)	Null models, adjusted for age, sex and an individual protein abundance.
# b)	Full models, that will include CVD risk factors (such as age, sex, systolic blood pressure,
#       total cholesterol, HDL-cholesterol, LDL-cholesterol, pack years of smoking, family history of CVD,
#       rheumatoid arthritis, diabetes, BMI, years of education, SIMD score) and individual protein abundance.
import os
import pyreadr
import pandas as pd
import numpy as np
from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test
from lib.pandas_ext import two_dfs_merge
import lib.stats as stats


flag = "hosp"  # hosp_gp, hosp, hosp_gp_cons
run = "agesex"
interesting_events = ["myocardial_infarction", "isch_stroke", "hf", "chd_nos",
                      "tia", "composite_CVD", "CVD_death"]
event = interesting_events[0]

proteins = pyreadr.read_r("data/phenotypes/GS_ProteinGroups_RankTransformed_23Aug2023.rds")[None]
proteins.set_index("id", inplace=True)
path = f"results/incremental_models/{run}/{flag}"
cox_path = f"results/cox/{flag}/{flag}_{event}.csv"

if not os.path.exists(path):
    os.makedirs(path)
    print(f"Path: {path} created!")

# for event in interesting_events:

cox = pd.read_csv(cox_path)
# this dataset requires more risk factors (family history of CVD left) - as mentioned above
cox = cox[['id', 'age', 'sex', 'avg_sys', 'Total_cholesterol',
           'HDL_cholesterol', 'pack_years', 'rheum_arthritis_Y', 'diabetes_Y',
           'bmi', 'years', 'rank', 'event', 'tte']]  # 17529
cox.set_index("id", inplace=True)
cox['sex'] = cox['sex'].replace({'M': 1, 'F': 0})

# %%
# so much work just to find out I need a merged dataframe xD
cox, proteins = two_dfs_merge(cox, proteins)

sc = StandardScaler(with_mean=True)
sc.fit(proteins)
sc.scale_ = np.std(proteins, axis=0, ddof=1).to_list()
sc.transform(proteins)

hazard_ratios = pd.DataFrame()
concordance = pd.DataFrame()

protein_name = "Total_cholesterol"

cph = CoxPHFitter()
cph.fit(cox, duration_col='tte', event_col='event', formula=f"age+sex")
test = proportional_hazard_test(cph, cox, time_transform="rank")
# R uses the default `km`, we use `rank`, as this performs well versus other transforms.
 ## age does not meet PH assumptions


# %%


# %%

#
# path = paste0('Lab/Proteins/individual_assoc/runs/deaths_included/', run, '/', flag)
# dir.create(path)
# # event = interesting_events[1]
# for (event in interesting_events) {

#
# null = coxph(Surv(merged$tte, merged$event) ~
# scale(transform(merged$score2))
# )
#
# i = 1
#
# output = paste0(path, '/', event)
# dir.create(output)
# mod = coxph(Surv(merged$tte, merged$event) ~
# scale(transform(merged$score2)) +
# scale(merged[, start]))
#
# outfile = paste0(output, '/', "model_summary.txt")
# text = sprintf("n=%d, number of events=%d", summary(mod)$n, summary(mod)$nevent)
# cat(text, file=outfile)
#
# for (prot in merged[, start:end]) {
#
#     protein = colnames(merged)[start-1+i]
#
# mod = coxph(Surv(merged$tte, merged$event) ~
# scale(transform(merged$score2)) +
# scale(prot))
#
#
# summary = serialize_coxph_summary(mod, protein)
# hazard_ratios = rbind(hazard_ratios, summary)
#
# cidx = concordance(null, mod)
#
# c_null = cidx$concordance["null"]
# c_mod = cidx$concordance["mod"]
# diff = c_mod - c_null
#
# concordance = rbind(concordance, data.frame(protein, c_null, c_mod, diff))
#
# i = i+1
# }
#
# f1 = paste0(output, '/HRs_score2_proteins_', event, '.csv')
# f2 = paste0(output, '/concordance_score2_proteins_', event, '.csv')
#
# write.csv(hazard_ratios, f1, row.names = F)
# write.csv(concordance, f2, row.names = F)
#
# hazard_ratios = subset(hazard_ratios, p < 0.05)
#
# f1 = paste0(output, '/HRs_score2_proteins_', event, '_significant.csv')
#
# write.csv(hazard_ratios, f1, row.names = F)
# }


#
# extract_coxme_table < - function(mod)
# {
#     beta < - mod$coefficients  # $fixed is not needed
# hr < -exp(beta)
# nvar < - length(beta)
# nfrail < - nrow(mod$var) - nvar
# se < - sqrt(diag(mod$var)[nfrail + 1: nvar])
# z < - round(beta / se, 2)
# p < - signif(1 - pchisq((beta / se) ^ 2, 1), 2)
# lci = exp(beta - 1 * 1.96 * se)
# uci = exp(beta + 1 * 1.96 * se)
# table = data.frame(cbind(beta, hr, se, z, lci, uci, p))
# return (table)
# }
#
# # Rank Based Inverse Normalisation of the data
# transform < - function(x)
# {
# transformed < - qnorm((rank(x, na.last="keep") - 0.5) / sum(! is.na(x)))
# return (transformed)
# }
#
# serialize_coxph_summary < - function(mod, protein)
# {
# coefs = extract_coxme_table(mod)
# hr < -coefs[2, 2]
# p < -coefs[2, 7]
# lci < -coefs[2, 5]
# uci < -coefs[2, 6]
# assumptions = cox.zph(mod)
# local = assumptions$table["scale(prot)", 'p']
# global=assumptions$table["GLOBAL", 'p']
# ret = data.frame(hr, p, lci, uci, local,
# global, protein)
# return (ret)
# }
#