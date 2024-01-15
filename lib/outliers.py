import numpy as np
import pandas as pd
from lib.stats import scale


def outlier_id(x, cut=4):
    xx = scale(x)
    return abs(xx) >= cut


def outlier_trim(x, cut=4):
    out_id = outlier_id(x, cut)
    tmp = x.copy()
    tmp[out_id] = np.NaN

    return tmp

