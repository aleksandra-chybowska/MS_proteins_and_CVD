from unittest import TestCase
import pandas as pd
import pyreadr
from lib.pandas_ext import two_dfs_merge


class TestPandasExt(TestCase):
    def test_two_dfs_merge(self):

        proteins = pyreadr.read_r("../data/phenotypes/GS_ProteinGroups_RankTransformed_23Aug2023.rds")[None]
        cox = pd.read_csv("../results/cox/hosp/hosp_myocardial_infarction.csv")

        merged = pd.merge(proteins, cox, left_index=True, right_index=True)
        cox, proteins = two_dfs_merge(cox, proteins)

        self.assertTrue(cox.index.equals(proteins.index))
        self.assertTrue(merged.index.equals(proteins.index))
