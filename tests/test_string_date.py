import unittest

import numpy as np

from lib.string_date import date_diff


class StringDateTest(unittest.TestCase):
    def test_date_diff(self):
        date1 = "200703"
        date2 = "201305"
        tte = date_diff(date1, date2)

        self.assertAlmostEqual(tte, 6.167, places=3)

        date3 = "200708"
        tte = date_diff(date3, date2)
        self.assertAlmostEqual(tte, 5.750)

        date4 = np.NaN
        tte = date_diff(date1, date4)
        np.testing.assert_equal(tte, np.NaN)

