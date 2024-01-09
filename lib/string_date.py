# Module to deal with yyyymm dates from GS
import math
import numpy as np


def date_to_year_month(date):
    d = str(date)
    return d[0:4], d[4:6]


def year_month_to_date(year, month):
    return f"{int(year)}{int(month):02d}"


def date_diff(baseline_date, event_date):
    """
    :param baseline_date: yyyymm date - in Cox analysis, this is the baseline appt
    :param event_date: yyyymm date - in Cox analysis, this is the event
    :return: difference between the two dates (in years)
    """
    if math.isnan(float(baseline_date)) or math.isnan(float(event_date)):
        return np.NaN

    appt_year = int(str(baseline_date)[0:4])
    appt_month = int(str(baseline_date)[4:6])

    dis_year = int(str(event_date)[0:4])
    dis_month = int(str(event_date)[4:6])

    tte = (dis_year - appt_year) + ((dis_month - appt_month) / 12)

    return tte

