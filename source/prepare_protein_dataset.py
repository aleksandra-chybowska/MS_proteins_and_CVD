import pyreadr
from sklearn.preprocessing import StandardScaler

from helpers.womans_health import get_on_pill
from lib.parquet_helper import write_parquet
from lib.stats import *

pheno = pyreadr.read_r("data/phenotypes/GS_phenos_internal_with_DST_28Nov2023_REM.rds")[None]
proteins = pyreadr.read_r("data/phenotypes/GS_ProteinGroups_RankTransformed_23Aug2023.rds")[None]
annots = pd.read_csv("data/annotations/short_annots.csv")
pill = pd.read_table("data/disease/womans_phenotypes/GS_womens_phenotypes_v2v5combined.txt")


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

# scaling
cols = proteins.columns[1:]
females[cols] = scale(females[cols])

# %%
path = "results/incremental_models/plot_data/females_and_proteins_scaled.parquet"
write_parquet(females, path)
print(f"Dataset prepared: {path}")


# quick test - it was ok

# sc = StandardScaler(with_mean=True)
# sc.fit(females[cols])
# sc.scale_ = np.std(females[cols], axis=0, ddof=1).to_list()
# sklearn_object = sc.transform(females[cols])
# sklearn_object = pd.DataFrame(sklearn_object)
# np.allclose(females[cols[0]].values, sklearn_object.iloc[:, 0].values)
