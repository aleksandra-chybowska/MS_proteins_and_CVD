# %%
import pyreadr
import pandas as pd
import os
from lib.parquet_helper import read_parquet
from plotnine import ggplot, aes, labs, geom_boxplot

pheno = pyreadr.read_r("data/phenotypes/GS_phenos_internal_with_DST_28Nov2023_REM.rds")[None]
proteins = read_parquet("data/transformed_input/GS_ProteinGroups_RankTransformed_23Aug2023_scaled.parquet")

general = pd.merge(pheno, proteins, on="id")

# %%
plot = (ggplot(general, aes(x="sex", y="P01019")) +
        labs(title="Angiotensinogen by sex - no age limit", y="Angiotensinogen") +
        geom_boxplot()).draw(False)
plot.show()
plot.savefig("plots/self_reported_pill_full_ds.pdf")

# %%
plot = (ggplot(general.loc[(general.age < 50) & (general.age > 18)], aes(x="sex", y="P01019")) +
        labs(title="Angiotensinogen by sex - 18-50yrs", y="Angiotensinogen") +
        geom_boxplot()).draw(False)
plot.show()
plot.savefig("plots/self_reported_pill_18age50.pdf")

# %%
females = read_parquet("results/incremental_models/plot_data/females_and_proteins.parquet")
females["on_pill"] = pd.Categorical(females["on_pill"])
plot = (ggplot(females, aes(x="on_pill", y="P01019")) +
        labs(title="Angiotensinogen by pill - no age limit", y="Angiotensinogen") +
        geom_boxplot()).draw(False)
plot.show()
plot.savefig("plots/angiotensinogen_by_pill_full.pdf")

# %%
plot = (ggplot(females.loc[(females.age < 50) & (females.age > 18)], aes(x="on_pill", y="P01019")) +
        labs(title="Angiotensinogen by pill - 18-50yrs", y="Angiotensinogen") +
        geom_boxplot()).draw(False)
plot.show()
plot.savefig("plots/angiotensinogen_by_pill_18age50.pdf")

# %%
# compare females on pill and males in the same age
