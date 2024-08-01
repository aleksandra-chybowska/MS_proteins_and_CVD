#%%
import pandas as pd
from upsetplot import UpSet
from collections import Counter
import matplotlib.pyplot as plt

flag = "hosp"
type = "40-69"
run = "agesex"

interesting_events = ["myocardial_infarction", "isch_stroke", "hf", "chd_nos",
                      "tia", "death", "CVD_death", "composite_CVD"]

all_data = pd.DataFrame()
for event in interesting_events:
    csv = pd.read_csv(f"results/cox/{type}/cox_{flag}_{event}_prepped.csv")
    csv["event_name"] = event
    all_data = pd.concat([all_data, csv])

all_data = all_data.query("event == 1")
all_data = all_data[["id", "event_name", "event"]]

all_data_pivoted = all_data.pivot(index=['id'], columns="event_name", values="event")
all_data_pivoted = all_data_pivoted.fillna(0)
all_data_pivoted.to_csv("plots/upset_as_table.csv")

# Create a binary matrix for the combinations
data_binary = all_data_pivoted.astype(bool)

# Count the unique rows (combinations of events)
combination_counts = Counter([tuple(row) for row in data_binary.values])

# Convert the Counter to a Pandas Series
combinations = pd.Series(combination_counts).sort_index()

# Convert the index to MultiIndex
combinations.index = pd.MultiIndex.from_tuples(combinations.index, names=all_data_pivoted.columns)
combinations.to_csv("plots/upset_as_table.csv")

# Print the series to inspect
print(combinations)
upset = UpSet(combinations, facecolor="darkblue", show_counts="{:,}", sort_by="cardinality")
upset.plot()
#%%
plt.suptitle('UpSet Plot of Event Overlaps')
plt.xlabel('Event Combinations')
plt.ylabel('Number of Individuals')

upset.plot()
plt.savefig(f"plots/upset.png")
