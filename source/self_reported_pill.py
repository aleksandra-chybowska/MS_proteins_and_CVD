# %%
import pyreadr
import pandas as pd
import numpy as np
import statsmodels.api as sm
from lib.lm import results_summary_to_dataframe
from helpers.womans_health import get_on_pill
from plotnine import ggplot, aes, geom_boxplot, labs, ggsave

pheno = pyreadr.read_r("data/phenotypes/GS_phenos_internal_with_DST_28Nov2023_REM.rds")[None]
proteins = pyreadr.read_r("data/phenotypes/GS_ProteinGroups_RankTransformed_23Aug2023.rds")[None]
annots = pd.read_csv("data/annotations/short_annots.csv")
pill = pd.read_table("data/disease/womans_phenotypes/GS_womens_phenotypes_v2v5combined.txt")

# %%
general = pd.merge(pheno, proteins, on="id")
plot = (
   ggplot(general)
   + aes(x="sex", y="P01019")
   + labs(title="Sex differences in Angiotensinogen levels - no age limit", y="Angiotensinogen")
   + geom_boxplot()
)

ggsave(plot, filename="plots/self_reported_pill_full_ds.pdf")

# typical age of menopause between 45 and 55. In many Western countries it is 51 yrs.
plot = (
   ggplot(general.loc[(general.age < 50) & (general.age > 18)])
   + aes(x="sex", y="P01019")
   + labs(title="Sex differences in Angiotensinogen levels - age between 18 and 50", y="Angiotensinogen")
   + geom_boxplot()
)

ggsave(plot, filename="plots/self_reported_pill_18age50.pdf")

# %%
females = pd.merge(pheno, pill, on="id")
females = females[females.sex == "F"]
females = females[~np.isnan(females.taken_cont)]
females.loc[:, "taken_cont"] = females["taken_cont"].apply(lambda x: 1 if x == 1 else 0)

cols = ["id", "age", "sex", "taken_cont", "age_started_cont", "years_taking_cont"]
females = females[cols]
females["on_pill"] = females.apply(
    lambda row: get_on_pill(row.taken_cont, row.age, row.age_started_cont, row.years_taking_cont), axis=1)
females = pd.merge(females, proteins, on="id")

pd.crosstab(columns="count", index=females["on_pill"])  # 2082 females that could be on pill

# quick checks
# angiotensinogen plotted for males and females
ggplot(females) + aes(x="sex", y="P01019") + geom_boxplot()
# angiotensinogen plotted for females on pill and not on pill

# pca of angiotensinogen



cols.append("on_pill")
females_pheno = females[cols]
females_proteins = females[females.columns.difference(cols)]

# %%
# standardise proteins
#

#%%
df = pd.DataFrame()
for protein in females_proteins.columns:

    # define response variable
    y = females_proteins[[protein]]
    x = sm.add_constant(females_pheno[["on_pill", "age"]])
    feature = "on_pill"
    # fit linear regression model
    model = sm.OLS(y, x).fit()
    summary = results_summary_to_dataframe(model)
    new_data = pd.DataFrame({"Protein": [protein],
                             "Coef": [summary.loc[feature, 'coeff']],
                             "P": [summary.loc[feature, 'pvals']]})
    df = pd.concat([df, new_data], ignore_index=True)

df = pd.merge(df, annots, left_on="Protein", right_on="id", right_index=False).drop(columns='id')
# %%


df = df.sort_values(by='P', ascending=True)
df.loc[df.P<0.05]
# divide by the number of experiments - bonferroni correction



df.to_csv("results/proteins_that_vary_with_pill.csv", index=False)
