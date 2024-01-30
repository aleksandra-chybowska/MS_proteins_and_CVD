# %%
import matplotlib.pyplot as plt
import pandas as pd

analysis_type = "plotly_express"
df = pd.read_csv("results/incremental_parallel/hosp/agesex/merged_results.csv")
events = df["event"].unique()
formula_basic = "age+sex+protein"
formula_full = ("age+sex+protein+avg_sys+Total_cholesterol+HDL_cholesterol+"
                "pack_years+rheum_arthritis_Y+diabetes_Y+years+rank+on_pill")
bonf_correction = 0.05/439
threshold = bonf_correction
# df = df[df["id"] == "P01019"]
plotting = pd.DataFrame()
# %%
for event in events:
    tmp = df[df["event"] == event]
    basic = tmp[tmp["formula"] == formula_basic]
    basic = basic[["id", "name", "hr", "lci", "uci", "p"]]
    basic.rename(columns={"hr": "hr_basic", "lci": "lci_basic", "uci": "uci_basic", "p": "p_basic"},
                 inplace=True)
    basic["err_basic"] = (basic["uci_basic"] - basic["lci_basic"]) / 2
    significant_basic = basic[basic["p_basic"] < threshold]
    names = set(significant_basic["name"])

    full = tmp[tmp["formula"] == formula_full]
    full = full[["id", "hr", "lci", "uci", "event", "name", "p"]]
    full.rename(columns={"hr": "hr_full", "lci": "lci_full", "uci": "uci_full", "p": "p_full"},
                inplace=True)
    full["err_full"] = (full["uci_full"] - full["lci_full"]) / 2
    significant_full = full[full["p_full"] < threshold]
    names.update(set(significant_full["name"]))

    tmp = pd.merge(basic, full, on="id", how="inner")
    tmp = tmp.query("name_x in @names").copy()
    tmp["both_significant"] = tmp.apply(lambda x:
                                        "yes" if x["p_basic"] < threshold and x["p_full"] < threshold else "no", axis=1)
    plotting = pd.concat([plotting, tmp], axis="rows")

# %%
if analysis_type == "plotly_express":
    import plotly.express as px
    import plotly.io as pio

    pio.renderers.default = "browser"
    fig = px.scatter(plotting, x="hr_basic", y="hr_full", opacity=0.7, color="both_significant",
                     error_x="err_basic", error_y="err_full", hover_name="name_x",
                     facet_col="event", facet_col_wrap=2,
                     labels={
                         "hr_basic": "HR basic (95% CI)",
                         "hr_full": "HR full (95% CI)",
                         "event": "Event",
                         "both_significant": "Both HRs significant (Bonferroni)"
                     })
    fig.for_each_annotation(lambda a: a.update(text=a.text.split("=")[-1]))
    fig.show()
elif analysis_type == "matplotlib":
    plt.errorbar(plotting["hr_basic"], plotting["hr_full"],
                 xerr=plotting["err_basic"],
                 yerr=plotting["err_full"], fmt='o', alpha=0.6)
    plt.xlabel("HR basic")
    plt.ylabel("HR full")
    plt.title("")
    plt.show()
else:
    from plotly.subplots import make_subplots
    import plotly.graph_objects as go

    fig = make_subplots(
        rows=2, cols=3,
        subplot_titles=("Plot 1", "Plot 2", "Plot 3", "Plot 4", "Plot 5", "Plot 6"))

    fig.add_trace(go.Scatter(x=plotting["hr_basic"], y="hr_full", opacity=0.7,
                             error_x="err_basic", error_y="err_full", hover_name="name_x"),
                  row=1, col=1)

    fig.add_trace(go.Scatter(x=[20, 30, 40], y=[50, 60, 70]),
                  row=1, col=2)

    fig.add_trace(go.Scatter(x=[300, 400, 500], y=[600, 700, 800]),
                  row=1, col=3)

    fig.add_trace(go.Scatter(x=[4000, 5000, 6000], y=[7000, 8000, 9000]),
                  row=2, col=1)

    fig.add_trace(go.Scatter(x=[300, 400, 500], y=[600, 700, 800]),
                  row=2, col=2)

    fig.add_trace(go.Scatter(x=[4000, 5000, 6000], y=[7000, 8000, 9000]),
                  row=2, col=3)

    fig.update_layout(height=900, width=1600,
                      title_text="Multiple Subplots with Titles")

    fig.show()
