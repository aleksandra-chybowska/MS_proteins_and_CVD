import pandas as pd
import statsmodels.api as sm


def results_summary_to_dataframe(results):
    """
    take the result of an statsmodel results table and transforms it into a dataframe
    :param results:
    :return: a dataframe with results of lm (R lm summary alike)
    """
    pvals = results.pvalues
    coef = results.params
    conf_lower = results.conf_int()[0]
    conf_higher = results.conf_int()[1]

    results_df = pd.DataFrame({"pvals": pvals,
                               "coef": coef,
                               "conf_lower": conf_lower,
                               "conf_higher": conf_higher
                                })

    # Reordering...
    results_df = results_df[["coef","pvals","conf_lower","conf_higher"]]
    return results_df

