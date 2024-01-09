import unittest
import helpers.womans_health as wh


class MyTestCase(unittest.TestCase):
    def test_translate_categories_age_started(self):
        self.assertEqual(wh.translate_categories_age_started(1)[0], 0)
        self.assertEqual(wh.translate_categories_age_started(2)[0], 20)
        self.assertEqual(wh.translate_categories_age_started(2)[1], 24)
        self.assertEqual(wh.translate_categories_age_started(3)[0], 25)
        self.assertEqual(wh.translate_categories_age_started(3)[1], 29)
        self.assertEqual(wh.translate_categories_age_started(4)[0], 30)
        self.assertEqual(wh.translate_categories_age_started(4)[1], 34)

    def test_get_on_pill(self):
        taken_cont = 0
        age = 26
        age_started = 1  # true age: 18, category == 2
        duration = 4  # true duration: 8 yrs, category == 4

        self.assertEqual(wh.get_on_pill(taken_cont, age, age_started, duration), 0)

        taken_cont = 1
        self.assertEqual(wh.get_on_pill(taken_cont, age, age_started, duration), 1)

        age = 50
        age_started = 2  # 20 - 24
        duration = 3  # between 3 and 4 years
        self.assertEqual(wh.get_on_pill(taken_cont, age, age_started, duration), 0)
