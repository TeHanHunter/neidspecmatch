import pandas as pd
import matplotlib.pyplot as plt

def plot_vsini_distribution_neid(csv_file, orders=None):
    """
    Plot vsini distribution for stars in the CSV file.

    Parameters:
    -----------
    csv_file : str
        Path to the CSV file containing vsini data
    orders : list, optional
        List of orders to plot. If None, all orders will be plotted.
    """
    df = pd.read_csv(csv_file, sep=",")
    df["TIC_ID"] = df["TIC"].str.replace("TIC ", "").astype(str)
    df["Observation_ID"] = df["Folder"] + "_" + df["TIC_ID"]

    all_vsini_cols = ["vsini_55", "vsini_101", "vsini_102", "vsini_103"]
    vsini_cols = orders if orders is not None else all_vsini_cols

    df_long = df.melt(
        id_vars=["TIC_ID", "Observation_ID", "Folder"],
        value_vars=vsini_cols,
        var_name="Order",
        value_name="vsini"
    )

    medians = df_long.groupby("Observation_ID")["vsini"].median().reset_index()
    df_long = df_long.merge(medians.rename(columns={"vsini": "vsini_median"}), on="Observation_ID")

    tic_medians = df_long.groupby("TIC_ID")["vsini"].median().reset_index()
    sorted_tics = tic_medians.sort_values("vsini")["TIC_ID"].tolist()
    df_long["TIC_ID"] = pd.Categorical(df_long["TIC_ID"], categories=sorted_tics, ordered=True)

    order_colors = {
        "vsini_55": "tab:blue",
        "vsini_101": "tab:orange",
        "vsini_102": "tab:green",
        "vsini_103": "tab:red"
    }

    fig, ax = plt.subplots(figsize=(12, 6))
    for spine in ax.spines.values():
        spine.set_visible(True)

    orders_in_legend = set()

    for i, tic_id in enumerate(sorted_tics):
        tic_data = df_long[df_long["TIC_ID"] == tic_id]
        obs_ids = sorted(tic_data["Observation_ID"].unique())
        for j, obs_id in enumerate(obs_ids):
            obs_data = tic_data[tic_data["Observation_ID"] == obs_id]
            marker = "s" if j == 0 else "o"

            # Sort by order for consistent line plotting
            obs_data = obs_data.sort_values("Order")

            # Connect the vsini values of the same observation
            ax.plot([i]*len(obs_data), obs_data["vsini"],
                    color="gray", linestyle="-", zorder=1)

            for _, row in obs_data.iterrows():
                label = row["Order"] if row["Order"] not in orders_in_legend else ""
                if label:
                    orders_in_legend.add(row["Order"])

                # Hollow marker for the individual vsini
                ax.scatter(i, row["vsini"],
                           facecolors="none",
                           edgecolors=order_colors[row["Order"]],
                           marker=marker,
                           s=60,
                           label=label,
                           zorder=2)

            # Solid marker for the observation median
            ax.scatter(i, obs_data["vsini_median"].iloc[0],
                       color="black",
                       marker=marker,
                       s=60,
                       zorder=3)

    ax.set_xticks(range(len(sorted_tics)))
    ax.set_xticklabels(sorted_tics, rotation=90)
    ax.set_xlabel("TIC ID")
    ax.set_ylabel("vsini")
    ax.set_title("vsini Distribution by TIC ID")

    ax.legend(title="Order", loc="upper left")

    plt.tight_layout()
    plt.savefig("/Users/tehan/Documents/SURFSUP/NEID/NEIDSM/plots/vsini_distribution_by_tic.pdf",
                dpi=300, bbox_inches="tight")
    plt.show()


def plot_vsini_distribution_hpf(csv_file):
    """
    Plot vsini distribution ignoring empty rows.
    """

    # Read TSV (tab-separated) file
    df = pd.read_csv(csv_file, sep=",")
    df = df.dropna(how="all")
    df = df[df["TIC"].notna()]
    df = df[df["Folder"].notna()]
    vsini_cols = [col for col in df.columns if col.startswith("vsini_")]
    # Prepare TIC IDs
    df["TIC_ID"] = df["TIC"].str.replace("TIC_", "").astype(str)
    df["Observation_ID"] = df["Folder"] + "_" + df["TIC_ID"]

    # Reshape long
    df_long = df.melt(
        id_vars=["TIC_ID", "Observation_ID", "Folder"],
        value_vars=vsini_cols,
        var_name="Order",
        value_name="vsini"
    )

    # Drop rows where vsini is missing
    df_long = df_long.dropna(subset=["vsini"])

    # Compute observation medians
    medians = df_long.groupby("Observation_ID")["vsini"].median().reset_index()
    df_long = df_long.merge(medians.rename(columns={"vsini": "vsini_median"}), on="Observation_ID")

    # Compute TIC medians for sorting
    tic_medians = df_long.groupby("TIC_ID")["vsini"].median().reset_index()
    sorted_tics = tic_medians.sort_values("vsini")["TIC_ID"].tolist()
    df_long["TIC_ID"] = pd.Categorical(df_long["TIC_ID"], categories=sorted_tics, ordered=True)

    # Assign colors per order
    order_colors = dict(zip(sorted(vsini_cols), plt.cm.tab10.colors))

    fig, ax = plt.subplots(figsize=(12, 6))
    for spine in ax.spines.values():
        spine.set_visible(True)

    orders_in_legend = set()

    for i, tic_id in enumerate(sorted_tics):
        tic_data = df_long[df_long["TIC_ID"] == tic_id]
        obs_ids = sorted(tic_data["Observation_ID"].unique())

        for j, obs_id in enumerate(obs_ids):
            obs_data = tic_data[tic_data["Observation_ID"] == obs_id]
            marker_styles = ["s", "o", "X", "D", "^", "v", "p", "*", "h", "H", "<", ">"]
            marker = marker_styles[j % len(marker_styles)]

            # Connect lines within the same observation
            ax.plot(
                [i]*len(obs_data),
                obs_data["vsini"],
                color="gray",
                linestyle="-",
                zorder=1
            )

            for _, row in obs_data.iterrows():
                label = row["Order"] if row["Order"] not in orders_in_legend else ""
                if label:
                    orders_in_legend.add(row["Order"])

                ax.scatter(
                    i,
                    row["vsini"],
                    facecolors="none",
                    edgecolors=order_colors[row["Order"]],
                    marker=marker,
                    s=60,
                    label=label,
                    zorder=2
                )

            # Solid marker for median
            ax.scatter(
                i,
                obs_data["vsini_median"].iloc[0],
                color="black",
                marker=marker,
                s=40,
                zorder=3
            )

    ax.set_xticks(range(len(sorted_tics)))
    ax.set_xticklabels(sorted_tics, rotation=90)
    ax.set_xlabel("TIC ID")
    ax.set_ylabel("vsini")
    ax.set_title("vsini Distribution by TIC ID")

    ax.legend(title="Order", loc="upper left")

    plt.tight_layout()
    plt.savefig("/Users/tehan/Documents/SURFSUP/HPF/vsini_distribution_by_tic.pdf", dpi=300, bbox_inches="tight")
    plt.show()


if __name__ == '__main__':
    # plot_vsini_distribution_neid(
    #     "/Users/tehan/Documents/SURFSUP/NEID/NEIDSM/results_summary.csv",
    #     orders=["vsini_55", "vsini_101", "vsini_102", "vsini_103"]
    # )
    # plot_vsini_distribution_neid(
    #     "/Users/tehan/Documents/SURFSUP/NEID/NEIDSM/results_summary.csv",
    #     orders=["vsini_101", "vsini_102", "vsini_103"]
    # )
    plot_vsini_distribution_hpf(
        "/Users/tehan/Documents/SURFSUP/HPF/results_summary.csv",
    )
