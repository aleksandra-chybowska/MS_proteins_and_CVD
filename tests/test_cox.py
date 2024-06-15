import unittest
import numpy as np
from lib.cox import get_time_to_event


class MyTestCase(unittest.TestCase):
    def test_get_time_to_event(self):
        date_event = "201503"
        date_baseline = "200903"
        date_censor = "202204"
        date_death = "202108"

        # event happened - diff between baseline and event
        tte = get_time_to_event(date_baseline, date_event, date_censor, date_death)
        self.assertEqual(6, tte)

        date_event = np.NaN

        # in the event of death, this should return diff between baseline and death
        tte = get_time_to_event(date_baseline, date_event, date_censor, date_death)
        self.assertAlmostEqual(12.42, tte, places=2)

        date_death = np.NaN

        # no death and no event - return diff between baseline and the censor date
        tte = get_time_to_event(date_baseline, date_event, date_censor, date_death)
        self.assertAlmostEqual(13.08, tte, places=2)

        # weird cases
        # date of event is less than date of death

        date_death = "201503"
        date_baseline = "200903"
        date_censor = "202204"
        date_event = "202108"

        with self.assertRaises(Exception):
            get_time_to_event(date_baseline, date_event, date_censor, date_death)

        # This code has been removed as sometimes it is necessary to censor people with events
        # # date of event is greater than date of death
        # date_censor = "201503"
        # date_baseline = "200903"
        # date_event = "202204"
        # date_death = np.NaN
        #
        # with self.assertRaises(Exception):
        #     get_time_to_event(date_baseline, date_event, date_censor, date_death)

