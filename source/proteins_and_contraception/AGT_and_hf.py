# %%
import pandas as pd
from lifelines import CoxPHFitter
from lib.cox import extract_cox_coefs, summary_and_test

proteins = pd.read_csv(f"results/cox/40-69/proteins_hosp_all_events_scaled_8343.csv")
annots = pd.read_csv("data/annotations/short_annots.csv")

# %%
proteins = proteins[["id", "P01019"]]
df_hf = pd.read_csv("results/cox/40-69/cox_hosp_hf_prepped.csv") # 8343
cox = pd.merge(df_hf, proteins, left_on="id", right_on="id", how="inner") # 4874
cox['sex'] = cox['sex'].astype('category')
cox['on_pill'] = cox['on_pill'].astype('category')

# %%
pd.crosstab(columns="count", index=cox["on_pill"])  # 767 people that could be on pill

# Do I just run a full model (hf ~ age + sex + angiotensinogen + ... ) after excluding females with HCU?
# I’d add HCU as a binary var and see if it’s sig and if it alters the HR for agt.
# Could also look at a hcu:agt interaction.

cph = CoxPHFitter()
formula = "age + sex + P01019 + on_pill"

cph.fit(cox, duration_col='tte', event_col='event', formula=formula)
summary = cph.summary

# on pill is not significant in this model
#                     coef  exp(coef)  ...             p   -log2(p)
# covariate                            ...
# age             0.092490   1.096902  ...  4.808597e-13  40.919449
# sex[T.M]        0.512995   1.670286  ...  5.059020e-03   7.626926
# P01019         -0.397464   0.672022  ...  2.590755e-06  18.558196
# on_pill[T.1.0]  0.167660   1.182535  ...  6.294406e-01   0.667858

cph2 = CoxPHFitter()
formula2 = "age + sex + P01019"

cph2.fit(cox, duration_col='tte', event_col='event', formula=formula2)
summary2 = cph2.summary

#                coef  exp(coef)  se(coef)  ...         z             p   -log2(p)
# covariate                                 ...
# age        0.092316   1.096711  0.012779  ...  7.224188  5.041043e-13  40.851343
# sex[T.M]   0.486921   1.627298  0.173554  ...  2.805594  5.022393e-03   7.637409
# P01019    -0.396416   0.672726  0.084544  ... -4.688903  2.746740e-06  18.473848

cph3 = CoxPHFitter()
formula3 = "age + sex + P01019*on_pill"

cph3.fit(cox, duration_col='tte', event_col='event', formula=formula3)
summary3 = cph3.summary

#                            coef  exp(coef)  ...             p   -log2(p)
# covariate                                   ...
# age                    0.092421   1.096827  ...  4.870714e-13  40.900932
# sex[T.M]               0.510157   1.665553  ...  5.342118e-03   7.548372
# P01019                -0.406609   0.665905  ...  3.876728e-06  17.976729
# on_pill[T.1.0]         0.180296   1.197572  ...  6.047646e-01   0.725554
# P01019:on_pill[T.1.0]  0.115630   1.122581  ...  7.139983e-01   0.486007

cph4 = CoxPHFitter()
formula4 = "age + sex*P01019"

cph4.fit(cox, duration_col='tte', event_col='event', formula=formula4)
summary4 = cph4.summary

#                      coef  exp(coef)  ...             p   -log2(p)
# covariate                             ...
# age              0.092205   1.096589  ...  5.491026e-13  40.727990
# sex[T.M]         0.454906   1.576025  ...  1.502240e-02   6.056741
# P01019          -0.353296   0.702369  ...  4.593533e-03   7.766180
# sex[T.M]:P01019 -0.080448   0.922703  ...  6.356811e-01   0.653625
