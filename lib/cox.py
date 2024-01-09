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
        raise Exception("Something went wrong - we're censoring people with events.")

    if not math.isnan(float(date_event)):
        return tte

    if date_death is None:
        return tte_censor

    return min(tte_censor, tte_death)
