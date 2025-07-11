import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def plot_vsini_distribution(csv_file):
    df = pd.read_csv(csv_file, sep=",")  # Adjust separator if needed
    df["TIC_ID"] = df["TIC"].str.replace("TIC ", "").astype(str)

    vsini_cols = ["vsini_55", "vsini_101", "vsini_102", "vsini_103"]
    df_long = df.melt(id_vars="TIC_ID", value_vars=vsini_cols,
                      var_name="Order", value_name="vsini")

    medians = df_long.groupby("TIC_ID")["vsini"].median().reset_index()
    medians = medians.rename(columns={"vsini": "vsini_median"})
    df_long = df_long.merge(medians, on="TIC_ID")

    sorted_tics = medians.sort_values("vsini_median")["TIC_ID"].tolist()
    df_long["TIC_ID"] = pd.Categorical(df_long["TIC_ID"], categories=sorted_tics, ordered=True)

    order_colors = {
        "vsini_55": "tab:blue",
        "vsini_101": "tab:orange",
        "vsini_102": "tab:green",
        "vsini_103": "tab:red"
    }

    plt.style.use("seaborn-v0_8-darkgrid")
    fig, ax = plt.subplots(figsize=(10, 6))

    for i, tic in enumerate(sorted_tics):
        tic_data = df_long[df_long["TIC_ID"] == tic]
        for _, row in tic_data.iterrows():
            ax.scatter(i, row["vsini"], color=order_colors[row["Order"]], label=row["Order"] if i==0 else "", zorder=2)
        ax.plot([i, i], [tic_data["vsini"].min(), tic_data["vsini"].max()], color="gray", linestyle="--", zorder=1)
        ax.scatter(i, tic_data["vsini_median"].iloc[0], facecolors="none", edgecolors="black", s=80, zorder=3)

    ax.set_xticks(range(len(sorted_tics)))
    ax.set_xticklabels(sorted_tics, rotation=90)
    ax.set_xlabel("Star (TIC ID)")
    ax.set_ylabel("vsini")
    ax.set_title("vsini Distribution per Star with Median (Hollow Circle)")
    ax.legend(title="Order", loc="upper left")

    plt.tight_layout()
    plt.show()

if __name__ == '__main__':
    plot_vsini_distribution("/Users/tehan/Documents/SURFSUP/NEID/NEIDSM/results_summary.csv")