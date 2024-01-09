# "If so, how old were you when you first went on the contraceptive pill (...)?"
# "1 - <20, 2 - 20-24, 3 - 25-29, 4 - 30-34, 5 - 35-39, 6 - 40-44, 7 - 45-49, 8 - 50-54, 9 - 55-59, 10 - 60+"
def translate_categories_age_started(category):
    if category < 1 or category > 10:
        raise Exception("Category must be between 1 and 10")

    end = 20+5*(category-1)-1
    start = end-4

    if category == 1:
        return 0, 19
    elif 1 < category < 10:
        return start, end
    else:
        return 60, 100


# "For about how many years in total were you taking the contraceptive pill (...)"
# "1 - <1, 2 - 1-2, 3 - 3-4, 4 - 5-9, 5 - 10-14, 6 - 15-19, 7 - 20-24, 8 - 25+"
def translate_categories_duration(category):
    if category < 1 or category > 8:
        raise Exception("Category must be between 1 and 8")

    if category == 1:
        return 0, 1
    elif category == 2:
        return 1, 2
    elif category == 3:
        return 3, 4
    elif category == 4:
        return 5, 9
    elif category == 5:
        return 10, 14
    elif category == 6:
        return 15, 19
    elif category == 7:
        return 20, 24
    else:
        return 25, 80


# compare age to age started contraceptive therapy and for how many years they were taking
# age ~= (age started + for how long)
# 26 ~= 18 + 8
def get_on_pill(taken_cont, age, age_started, duration):
    # original encoding sets no as 2, I prefer to have zeros for no
    if taken_cont == 2 or taken_cont == 0:
        return 0

    window_min = translate_categories_age_started(age_started)[0] + translate_categories_duration(duration)[0]
    window_max = translate_categories_age_started(age_started)[1] + translate_categories_duration(duration)[1]

    return int(window_min < age < window_max)
