import pandas as pd
from lifelines.statistics import proportional_hazard_test
from matplotlib import pyplot as plt
import warnings
from lib.string_date import date_diff
import math


def get_time_to_event(date_baseline, date_event, date_censor, date_death=None):
    """
    Calculates the time to event for cox analysis
    :param date_baseline: recruitment date
    :param date_event: date of event
    :param date_censor: censor date
    :param date_death: date of death (if deaths are censored)
    :return:
    """

    # check for possible errors
    tte = date_diff(date_baseline, date_event)
    tte_death = date_diff(date_baseline, date_death)
    tte_censor = date_diff(date_baseline, date_censor)

    if math.isnan(float(tte_censor)):
        raise Exception("Time to censoring is NaN.")

    if date_death is not None and tte_death < tte:
        raise Exception("Event occurred after patient's death.")

    if tte_censor < tte:
        # Maybe a warning would work better here?
        print("Warning - we're censoring people with events.")
        return min(tte_censor, tte_death)

    if not math.isnan(float(date_event)):
        return tte

    if date_death is None:
        return tte_censor

    return min(tte_censor, tte_death)


def summary_and_test(cph, feature, df):
    """
    Extracts the most important results from lifeliness cox model summary
    :param cph: results of a lifeliness cox model
    :param feature: main parameter of the report (age? sex? protein?)
    :param df: dataframe containing variables to prepare cph model
    :return: Series with cox summary + cph test
    """

    cox_summary = cph.summary.loc[feature, :]
    variables = ['exp(coef)', 'p', 'exp(coef) lower 95%', 'exp(coef) upper 95%']
    new_names = ['hr', 'p', 'lci', 'uci']

    cox_summary = cox_summary[variables]
    cox_summary.index = new_names

    test = proportional_hazard_test(cph, df, time_transform="km")
    summary_test = test.summary.loc[feature, :]

    cox_summary["test_p"] = summary_test["p"]
    first = pd.Series({"feature": feature})
    ret = pd.concat([first, cox_summary])
    return ret


def extract_cox_coefs(cph, feature):
    """
    Extracts the most important results from lifeliness cox model summary for a variable
    :param cph: results of a lifeliness cox model
    :param feature: main parameter of the report (age? sex? protein?)
    :return: Series with cox summary
    """
    cox_summary = cph.summary.loc[feature, :]
    variables = ['coef', 'exp(coef)', 'se(coef)', 'exp(coef) lower 95%', 'exp(coef) upper 95%', 'z', 'p']
    new_names = ['beta', 'hr', 'se', 'lci', 'uci', 'z', 'p']

    cox_summary = cox_summary[variables]
    cox_summary.index = new_names

    first = pd.Series({"feature": feature})
    ret = pd.concat([first, cox_summary])

    return ret


def plot_partial_effects(cph, covariates, range, legend, save=False):
    """
    Plots partial effects of covariates on cox model.
    :param cph: results of a lifeliness cox model
    :param covariates: array containing covariates (sex? protein?). Example: ["sex", "protein"]
    :param range: array of arrays defining values of covariates. Example: [["M", -3], ["M", 3], ["F", -3], ["F", 3]]
    :param legend: how do we label these covariates + ranges on the plot
    :param save: if true, plots will be saved to "plots" dir
    :return: matplotlib plot
    """
    cph.plot_partial_effects_on_outcome(covariates=covariates,
                                        values=range, cmap='coolwarm', plot_baseline=False)
    plt.ylabel("HF-free survival")
    plt.xlabel("Follow up (years)")
    plt.legend(legend, loc="lower left")
    if save:
        plt.savefig(f"plots/partial_effects_{covariates[0]}_{covariates[1]}.png", dpi=600)
    plt.show()
