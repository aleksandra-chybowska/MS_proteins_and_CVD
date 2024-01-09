# %%
from sksurv.datasets import load_veterans_lung_cancer
from lifelines import CoxPHFitter
import pandas as pd
from sksurv.preprocessing import OneHotEncoder

data_x, data_y = load_veterans_lung_cancer()
data_y = pd.DataFrame(data_y)
df = pd.merge(data_x, data_y, left_index=True, right_index=True)

# %%
encoder = OneHotEncoder()
encoder.fit(df)
data_transformed = encoder.transform(df)

# %%
cph = CoxPHFitter()
cph.fit(data_transformed, duration_col='Survival_in_days', event_col='Status')
cph.print_summary()

result = cph.summary.T
