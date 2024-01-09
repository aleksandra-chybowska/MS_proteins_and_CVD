# %%
import pyreadr
from lib.stats import summary, scale
from lib.parquet_helper import *

proteins = pyreadr.read_r("data/phenotypes/GS_ProteinGroups_RankTransformed_23Aug2023.rds")[None]
proteins.set_index("id", inplace=True)
proteins_scaled = scale(proteins)

proteins_scaled.reset_index(drop=False, inplace=True)
write_parquet(proteins_scaled, "data/transformed_input/GS_ProteinGroups_RankTransformed_23Aug2023_scaled.parquet")

# quick test - it was ok

# sc = StandardScaler(with_mean=True)
# sc.fit(proteins)
# sc.scale_ = np.std(proteins, axis=0, ddof=1).to_list()
# sklearn_object = sc.transform(proteins)
# sklearn_object = pd.DataFrame(sklearn_object)
# (proteins_scaled.iloc[:, 4].values == sklearn_object.iloc[:, 4].values).all()
