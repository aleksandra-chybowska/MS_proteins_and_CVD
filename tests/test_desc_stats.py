import unittest
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from lib.stats import summary, scale


class DescStatsTest(unittest.TestCase):
    def test_summary(self):
        d = {'col1': [1, np.nan, np.nan], 'col2': [3, 4, np.nan]}
        d = pd.DataFrame(d)
        ret = summary(d)
        self.assertEqual(ret.loc['missing', 'col1'], 2)
        self.assertEqual(ret.loc['missing', 'col2'], 1)

    def test_scale(self):
        d = {'col1': [-2.24973459,  0.86352426, -0.20886251, 0.56296221, -1.88926, -1.16295529],
             'col2': [-2.24973459,  0.86352426, -0.20886251, 0.56296221, -1.88926, -1.16295529]}

        d = pd.DataFrame(d)
        ret = scale(d['col1'])

        sc = StandardScaler(with_mean=True)
        sc.fit(d)
        sc.scale_ = np.std(d, axis=0, ddof=1).to_list()
        sklearn_object = sc.transform(d)
        sklearn_object = pd.DataFrame(sklearn_object)

        r_result = [-1.2171302, 1.1979166, 0.3660345, 0.9647618, -0.9374994, -0.3740834]

        self.assertTrue((sklearn_object.iloc[:, 0].values == ret.values).all())
        self.assertTrue(np.allclose(ret.values, r_result))

        ret_df = scale(d)
        self.assertTrue(np.allclose(ret_df.col1.values, r_result))




