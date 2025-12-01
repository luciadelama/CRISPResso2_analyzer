import numpy as np
import matplotlib.pyplot as plt

def _grouped_barplot(df, value_col, title, ylabel, treat_order=None):
    # Required columns: "Treatment", "Replicate", and the metric (value_col)
    required = {"Treatment", "Replicate", value_col}
    if df.empty or not required.issubset(df.columns):
        return None

    # Group by Treatment and Replicate, compute mean value (just in case of duplicates)
    agg = (
        df.groupby(["Treatment", "Replicate"], as_index=False)[value_col]
          .mean()
    )
    replicates = agg["Replicate"].unique().tolist()
    treatments = agg["Treatment"].unique().tolist()
    if treat_order:
        present = treatments
        treatments = [t for t in treat_order if t in present] + [t for t in present if t not in treat_order]
    
    # Make sure the plot is wide enough to fit many treatments
    n_treat = len(treatments)
    n_rep = len(replicates)

    x = np.arange(n_treat)
    width = 0.8 / max(1, n_rep)

    fig, ax = plt.subplots(figsize=(max(6, n_treat*0.6), 3.5))
    colors = [plt.cm.tab10(i / n_rep) for i in range(n_rep)]

    # Loop over replicates to plot bars
    for i, rep in enumerate(replicates):
        y = []
        for t in treatments:
            m = agg.loc[(agg["Treatment"] == t) & (agg["Replicate"] == rep), value_col]
            y.append(float(m.iloc[0]) if len(m) else np.nan)
        offset = (i - (n_rep - 1) / 2) * width
        ax.bar(x + offset, y, width, label=str(rep), color=colors[i % len(colors)])

    # Set labels and title
    ax.set_xticks(x)
    ax.set_xticklabels(treatments, rotation=45, ha="right")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    if n_rep > 1:
        ax.legend(title="Replicate", ncols=min(n_rep, 3), frameon=False)
    fig.tight_layout()
    return fig

def make_frameshift_plot(df, treat_order=None):
    return _grouped_barplot(df, value_col="Frameshift%", title="Frameshift% per Treatment", ylabel="Frameshift %", treat_order=treat_order)

def make_indels_plot(df, treat_order=None):
    return _grouped_barplot(df, value_col="Indel%", title="Indels% per Treatment", ylabel="Indels %", treat_order=treat_order)

def make_sensitivity_plot(df, treat_order=None):
    # Required columns: "Treatment", "Replicate", and the metric (value_col)
    required = {"Treatment", "Replicate", "Sensitivity"}
    if df.empty or not required.issubset(df.columns):
        return None
    
    # Group by Treatment and Replicate, compute mean value (just in case of duplicates)
    agg = (
        df.groupby(["Replicate", "Treatment"], as_index=False)["Sensitivity"].mean()
    )
    replicates = agg["Replicate"].unique().tolist()
    treatments = agg["Treatment"].unique().tolist()
    if treat_order:
        present = treatments
        treatments = [t for t in treat_order if t in present] + [t for t in present if t not in treat_order]

    # Make sure the plot is wide enough to fit many replicates
    n_rep = len(replicates)
    n_treat = len(treatments)

    x = np.arange(n_rep)
    width = 0.8 / max(1, n_treat)

    fig, ax = plt.subplots(figsize=(max(6, n_rep*0.6), 3.5))
    colors = [plt.cm.Set2(v) for v in np.linspace(0, 1, max(1, n_treat))]

    # Loop over treatments to plot bars
    for i, t in enumerate(treatments):
        y = []
        for r in replicates:
            m = agg.loc[(agg["Replicate"] == r) & (agg["Treatment"] == t), "Sensitivity"]
            y.append(float(m.iloc[0]) if len(m) else np.nan)
        offset = (i - (n_treat - 1) / 2) * width
        ax.bar(x + offset, y, width, label=str(t), color=colors[i % len(colors)])

    # Set labels and title
    ax.set_xticks(x)
    ax.set_xticklabels(replicates, rotation=45, ha="right")
    ax.set_ylabel("Sensitivity")
    ax.set_title("Sensitivity per Replicate")
    if n_treat > 1:
        ax.legend(title="Treatment", ncols=min(n_treat, 3), frameon=False)
    fig.tight_layout()
    return fig