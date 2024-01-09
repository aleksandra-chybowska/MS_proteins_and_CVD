from unittest import TestCase

import numpy as np
import pandas as pd

from lib.int import rank_to_normal, rank_int_transform


class RankBasedInverseNormalTransformationTest(TestCase):
    def test_rank_to_normal(self):
        result = rank_to_normal(np.array([1, 2, 3, 4]), 0.5, 4)

        assert np.allclose(result,
                           np.array([-1.15034938, -0.31863936,  0.31863936,  1.15034938]))

    def test_rank_int_transform(self):
        s = pd.Series([2, 1, 1, np.NaN, 4, 3], index=["a", "b", "c", "d", "e", "f"])
        result = rank_int_transform(s, c=1.0/2)

        expected = pd.Series(index=['a', 'b', 'c', 'd', 'e', 'f'],
                             data=[0.0, -0.841621, -0.841621, np.nan, 1.281552, 0.524401])

        assert np.allclose(result, expected, atol=1e-3, equal_nan=True)

