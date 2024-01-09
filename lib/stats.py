import pandas as pd
import numpy as np


def summary(val):
    df = val.copy()
    if not isinstance(df, pd.DataFrame):
        df = pd.DataFrame(df)

    des1 = df.describe(include="all")
    des2 = df.isnull().sum().to_frame(name="missing").T
    return pd.concat([des1, des2])


def scale(y, c=True, sc=True):
    x = y.copy()

    if c:
        x -= x.mean()
    if sc and c:
        x /= x.std()
    elif sc:
        x /= np.sqrt(x.pow(2).sum().div(x.count() - 1))
    return x