import seaborn as sns
import pandas as pd

# Set seaborn style
sns.set(style="whitegrid")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
import seaborn as sns
sns.set(rc={'font.family': 'serif', 'font.serif': 'DejaVu Serif', 'font.size': 10,
            'axes.edgecolor': '0.2', 'axes.labelcolor': '0.', 'xtick.color': '0.', 'ytick.color': '0.',
            'axes.facecolor': '0.95', 'grid.color': '0.8'})

def plot_crossval_feh_delta_feh(feh_true, d_feh, ax=None):
    """
    INPUT:
        feh_true - True Fe/H
        d_feh - delta Fe/H = Fe/H_specmatch - Fe/H_True
    """
    if ax is None:
        fig, ax = plt.subplots(dpi=200)
    ax.plot(feh_true, d_feh, marker='o', lw=0, color='black', markersize=4)

    p = np.polyfit(feh_true, d_feh, deg=1)
    xx = np.linspace(-0.5, 0.5, 200)
    yy = np.polyval(p, xx)
    ax.plot(xx, yy, color='crimson', label='Linear fit\n$p_1$={:0.4f}\n$p_2$={:0.4f}'.format(p[0], p[1]))
    ax.legend(loc='upper right', fontsize='small')
    ax.set_xlabel('[Fe/H]')
    ax.set_ylabel('$\Delta$[Fe/H] = $\mathrm{[Fe/H]}_{\mathrm{Recovered}}$ - $\mathrm{[Fe/H]}_{\mathrm{True}}$', fontsize=8)
    ax.minorticks_on()


def plot_crossvalidation_results_main(order, df_crossval, savefolder=None):
    """
    Main crossvalidation results plot with seaborn styling and proper layout adjustments
    """
    PW = 10.
    PH = 3  # Adjust height to prevent overlapping
    fig = plt.figure(figsize=(PW, PH), dpi=300)
    gs0 = GridSpec(1, 2)
    gs0.update(top=0.92, bottom=0.11, wspace=0.05, left=0.07, right=0.62)  # Adjust left-right spacing

    gs1 = GridSpec(1, 1)
    gs1.update(top=0.92, bottom=0.11, hspace=0.4, left=0.7, right=0.98)

    ax = plt.subplot(gs0[0, 0])
    bx = plt.subplot(gs0[0, 1])
    cx = plt.subplot(gs1[0, 0])

    # Plotting Teff vs FeH
    ax.plot(df_crossval['feh_true'], df_crossval['teff_true'], marker='o', lw=0, markeredgecolor='black', markersize=4)

    # Plotting Teff vs logg with horizontal lines for clarity
    bx.plot(df_crossval['logg_true'], df_crossval['teff_true'], marker='o', lw=0, markeredgecolor='black', markersize=4)
    # Drawing connection lines for Teff vs FeH
    for i in range(len(df_crossval)):
        _x = [df_crossval['feh_true'].values[i], df_crossval['feh'].values[i]]
        _y = [df_crossval['teff_true'].values[i], df_crossval['teff'].values[i]]
        ax.plot(_x, _y, color='crimson', lw=0.5, zorder=1)

        _x = [df_crossval['logg_true'].values[i], df_crossval['logg'].values[i]]
        if i == 0:
            bx.plot(_x, _y, color='crimson', lw=0.5, label='Crossvalidation Value', zorder=1)
        else:
            bx.plot(_x, _y, color='crimson', lw=0.5, zorder=1)

    # Setting axis limits and labels with increased padding
    ax.set_xlim(-0.6, 0.6)
    ax.set_xlabel('[Fe/H]', labelpad=5)
    ax.set_ylabel('$T_{\mathrm{eff}}$ [K]', labelpad=5)
    bx.set_xlabel('$\log (g)$', labelpad=5)
    ax.minorticks_on()
    # bx.minorticks_on()
    bx.set_yticklabels([])  # Remove y-ticks for the middle plot to avoid clutter
    bx.set_ylim(*ax.get_ylim())

    fig.subplots_adjust(wspace=0.15)

    # Title for the plot
    fig.suptitle('Library Performance (NEID Order {})'.format(order), y=0.98, fontsize=10)

    # Call the delta FeH plot function
    plot_crossval_feh_delta_feh(df_crossval.feh_true.values, df_crossval.d_feh.values, ax=cx)

    # Save or display the figure
    if savefolder:
        fig.savefig('{}/crossvalidation_o{}_main.pdf'.format(savefolder, order), bbox_inches='tight')
    else:
        plt.show()


if __name__ == '__main__':
    order=102
    df_crossval = pd.read_csv(
        f'/Users/tehan/PycharmProjects/neidspecmatch/library/20240822_specmatch_nir/crossval/o{order}_crossval/crossvalidation_results_o{order}.csv')
    # Running the function with the provided data
    plot_crossvalidation_results_main(order=order, df_crossval=df_crossval,
                                      savefolder=f'/Users/tehan/PycharmProjects/neidspecmatch/library/20240822_specmatch_nir/crossval/o{order}_crossval/')
