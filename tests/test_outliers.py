import unittest

import pandas as pd
import numpy as np
from lib.parquet_helper import read_parquet
from lib.outliers import outlier_trim
from lib.stats import summary


class TestOutliers(unittest.TestCase):
    def test_outliers(self):
        #%%
        r = read_parquet("tests/mock/test_cholesterol_R.parquet")
        cox = read_parquet("tests/mock/test_cholesterol.parquet")
        cholesterol = cox["Total_cholesterol"]
        outliers = outlier_trim(cholesterol)
        outliers = pd.DataFrame(outliers).reset_index(drop=True).dropna()  # ty chuju
        r = r.dropna()

        # r[np.isnan(r["V1"])].index
        # outliers[np.isnan(outliers["Total_cholesterol"])].index
        self.assertTrue(np.allclose(r["V1"].values, outliers["Total_cholesterol"].values))  # add assertion here

