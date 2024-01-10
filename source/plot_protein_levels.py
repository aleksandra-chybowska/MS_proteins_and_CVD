# %%
import pyreadr
import pandas as pd
import os
from lib.parquet_helper import read_parquet
from plotnine import ggplot, aes, labs, geom_boxplot
from lib.pca import *


pheno = pyreadr.read_r("data/phenotypes/GS_phenos_internal_with_DST_28Nov2023_REM.rds")[None]
proteins = read_parquet("data/transformed_input/GS_ProteinGroups_RankTransformed_23Aug2023_scaled.parquet")

# %%
# general - proteins by sex
cols = pheno.columns
general = pd.merge(pheno, proteins, on="id")  # properly scaled proteins, phew!
general_pheno = general[pheno.columns]
general_proteins = general[general.columns.difference(cols)]

pca = run_PCA(general_proteins, 2)
pca_df_general = pd.concat([general_pheno, pca], axis=1)

# %%

# proteins - females only
females = read_parquet("results/incremental_models/plot_data/females_and_proteins_scaled.parquet")
females["on_pill"] = pd.Categorical(females["on_pill"])
cols = ["id", "age", "sex", "taken_cont", "age_started_cont", "years_taking_cont", "on_pill"]
females_pheno = females[cols]
females_proteins = females[females.columns.difference(cols)]

pca = run_PCA(females_proteins, 2)
pca_df = pd.concat([females_pheno, pca], axis=1)


# %%
plot = (ggplot(general, aes(x="sex", y="P01019")) +
        labs(title="Angiotensinogen by sex - no age limit", y="Angiotensinogen") +
        geom_boxplot()).draw(False)
plot.show()
plot.savefig("plots/self_reported_pill_full_ds.pdf")

# %%
# again, this is not entirely right - scaling was not needed. Not a big issue though
plot = (ggplot(general.loc[(general.age < 50) & (general.age > 18)], aes(x="sex", y="P01019")) +
        labs(title="Angiotensinogen by sex - 18-50yrs", y="Angiotensinogen") +
        geom_boxplot()).draw(False)
plot.show()
plot.savefig("plots/self_reported_pill_18age50.pdf")

# %%
plot = (ggplot(females, aes(x="on_pill", y="P01019")) +
        labs(title="Angiotensinogen by pill - no age limit", y="Angiotensinogen") +
        geom_boxplot()).draw(False)
plot.show()
plot.savefig("plots/angiotensinogen_by_pill_full.pdf")

# %%
# This may work better on unscaled data - scaling was not needed

plot = (ggplot(females.loc[(females.age < 50) & (females.age > 18)], aes(x="on_pill", y="P01019")) +
        labs(title="Angiotensinogen by pill - 18-50yrs", y="Angiotensinogen") +
        geom_boxplot()).draw(False)
plot.show()
plot.savefig("plots/angiotensinogen_by_pill_18age50.pdf")

# %%
# what to colour by?

target_col_name = "on_pill"
targets = [1, 0]
colors = ['r', 'b']
filename = "plots/PCA_females_on_pill_two_components.png"

plot_PCA(target_col_name, targets, colors, pca_df, filename)

target_col_name = "sex"
targets = ['F', 'M']
colors = ['r', 'b']
filename = "plots/PCA_proteins_by_sex_two_components.png"
plot_PCA(target_col_name, targets, colors, pca_df_general, filename)

