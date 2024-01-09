import pandas as pd
import statsmodels.api as sm


def results_summary_to_dataframe(results):
    """
    take the result of an statsmodel results table and transforms it into a dataframe
    :param results:
    :return:
    """
    pvals = results.pvalues
    coeff = results.params
    conf_lower = results.conf_int()[0]
    conf_higher = results.conf_int()[1]

    results_df = pd.DataFrame({"pvals": pvals,
                               "coeff": coeff,
                               "conf_lower": conf_lower,
                               "conf_higher": conf_higher
                                })

    # Reordering...
    results_df = results_df[["coeff","pvals","conf_lower","conf_higher"]]
    return results_df

