#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Function by EDM1, adapted by Aleksandra Chybowska
from copy import copy

import numpy as np
import pandas as pd
import scipy.stats as ss
from scipy.stats import norm


def rank_int_transform(series, c=1.0/2, stochastic=False):
    """ Perform rank-based inverse normal transformation on pandas series.
        If stochastic is True ties are given rank randomly, otherwise ties will
        share the same value. NaN values are ignored.

        Args:
            series (pandas.Series):   Series of values to transform
            c (Optional[float]): Constand parameter (Bloms constant = 3.0/8, Bliss constant = 1.0/2)
            stochastic (Optional[bool]):  Whether to randomise rank of ties

        Returns:
            pandas.Series
    """

    # Check input
    assert (isinstance(series, pd.Series))
    assert (isinstance(c, float))
    assert (isinstance(stochastic, bool))

    # Set seed
    np.random.seed(123)

    # Take original series indexes
    orig_idx = series.index

    orig = copy(series)
    # Drop NaNs
    series = series.loc[~pd.isnull(series)]

    # Get ranks
    if stochastic is True:
        # Shuffle by index
        series = series.loc[np.random.permutation(series.index)]
        # Get rank, ties are determined by their position in the series (hence
        # why we randomised the series)
        rank = ss.rankdata(series, method="ordinal")
    else:
        # Get rank, ties are averaged
        rank = ss.rankdata(series, method="average")

    # Convert numpy array back to series
    rank = pd.Series(rank, index=series.index)

    # Convert rank to normal distribution
    transformed = rank.apply(rank_to_normal, c=c, n=len(rank))
    orig[~np.isnan(orig)] = transformed
    orig[np.isnan(orig)] = np.nan
    return orig


def rank_to_normal(rank, c, n):
    # Standard quantile function
    x = (rank - c) / (n - 2 * c + 1)
    return ss.norm.ppf(x)


def test():
    # Unit tests in tests folder
    # Test
    s = pd.Series([2, 1, 1, np.NaN, 4, 3], index=["a", "b", "c", "d", "e", "f"])
    res = rank_int_transform(s, c=1.0/2)
    print(res)

    return 0


if __name__ == '__main__':
    test()
