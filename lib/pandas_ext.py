import pandas as pd
import numpy as np


def two_dfs_merge(indexed_df1, indexed_df2):
    """
    Manual merge of two dataframes based on their index. Two dataframes are returned.
    After running this function, the returned dataframes contain identical indexes.
    :param indexed_df1: indexed df 1 (usually index is placed on the id column)
    :param indexed_df2: indexed df 2
    :return: tuple, indexed dfs
    """

    ids = np.intersect1d(indexed_df1.index, indexed_df2.index)
    indexed_df1 = indexed_df1.loc[ids, :]
    indexed_df2 = indexed_df2.loc[ids, :]
    return indexed_df1, indexed_df2.reindex(ids)

