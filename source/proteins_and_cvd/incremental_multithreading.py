import gc
import logging
import multiprocessing as mp
from datetime import time
from typing import List, Dict
from multiprocessing.pool import ThreadPool as Pool

import pandas as pd

data_frame = pd.DataFrame()
samples = [
    {'id': 1, 'param1': 'val1', 'param2': 'val2', 'param3': 'val3'},
    {'id': 2, 'param1': 'val1', 'param2': 'val2', 'param3': 'val3'},
    {'id': 3, 'param1': 'val1', 'param2': 'val2', 'param3': 'val3'},
    {'id': 4, 'param1': 'val1', 'param2': 'val2', 'param3': 'val3'}
]


def calculate(samples: List[Dict]) -> List:
    with Pool(processes=mp.cpu_count()) as mp_pool:
        start = time.time()
        results = mp_pool.map(_calculate_cox, samples)
        gc.collect()
        end = time.time()
        logging.info("calculate time = %s" % (end - start))

    return results


def _calculate_cox(sample: Dict) -> Dict:
    """
    https://lifelines.readthedocs.io/en/latest/Survival%20Regression.html
    https://scikit-survival.readthedocs.io/en/stable/api/generated/sksurv.linear_model.CoxPHSurvivalAnalysis.html
    """
    # cox = SuperDuperCoxClass(lambda=sample['param1'])
    # cox.fit(data_frame)
    #
    # outcome = cox.predict()
    # time_to_survive = outcome['costam']
    #
    # return time_to_survive
    return True


def main():
    outcomes = calculate()
    results = pd.DataFrame(outcomes)
    results.to_csv('uber_fajne_rezultaty.csv')


if __name__ == "__main__":
    main()