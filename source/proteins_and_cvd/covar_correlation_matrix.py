# %%
import numpy as np
import pyreadr
import pandas as pd
import matplotlib.pyplot as plt
from lib.pandas_ext import two_dfs_merge
from scipy.stats import spearmanr

# %%
# 17150
pheno = pd.read_csv("data/transformed_input/generic_pheno.csv") # based on deaths_old, there's 1122 deaths in this file
proteins = pyreadr.read_r("data/phenotypes/GS_ProteinGroups_RankTransformed_23Aug2023.rds")[None]
pheno.set_index("id", inplace=True)
proteins.set_index("id", inplace=True)

pheno, proteins = two_dfs_merge(pheno, proteins)  # both 13374 records, need to recalculate events

# %%
covars = ['age', 'sex', 'avg_sys', 'Total_cholesterol',
          'HDL_cholesterol', 'pack_years', 'rheum_arthritis_Y', 'diabetes_Y',
          'bmi', 'years', 'rank', 'on_pill']
pheno_covars = pheno[covars]

# Apply mapping to the 'Sex' column
pheno_covars.loc[:, 'sex'] = pheno_covars['sex'].map({'M': 1.0, 'F': 0.0})
# calculate correlation matrix between covariates
corr_matrix = pheno_covars.corr(method="spearman")
corr_matrix.to_csv("results/correlations/correlation_matrix.csv", index=True)

#test
spearman_corr, p_value = spearmanr(pheno_covars['Total_cholesterol'], pheno_covars['HDL_cholesterol'])
spearman_corr  # 0.19

# %%
plt.figure(figsize=(11, 11))
plt.imshow(corr_matrix, cmap='coolwarm', interpolation='nearest')
plt.colorbar()
plt.title('Correlation Matrix')
plt.xticks(range(len(corr_matrix.columns)), corr_matrix.columns, rotation='vertical')
plt.yticks(range(len(corr_matrix.columns)), corr_matrix.columns)
plt.savefig("plots/corr_matrix.png")
plt.show()

# %%

