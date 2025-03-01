import seaborn as sns
import pandas as pd
import pickle

# Set seaborn style
sns.set(style="whitegrid")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
import seaborn as sns

sns.set(rc={'font.family': 'serif', 'font.serif': 'DejaVu Serif', 'font.size': 10,
            'axes.edgecolor': '0.2', 'axes.labelcolor': '0.', 'xtick.color': '0.', 'ytick.color': '0.',
            'axes.facecolor': '0.95', 'grid.color': '0.8'})
plt.rcParams['font.family'] = 'serif'
plt.rcParams['mathtext.fontset'] = 'dejavuserif'  # Use Computer Modern (serif font)


def plot_crossval_feh_delta_feh(feh_true, d_feh, ax=None,
                                tick_fontsize=8, ticklabel_fontsize=10,
                                xlabel_fontsize=12, ylabel_fontsize=12,
                                legend_fontsize=10):
    """
    INPUT:
        feh_true - True Fe/H
        d_feh - delta Fe/H = Fe/H_specmatch - Fe/H_True
    """
    if ax is None:
        fig, ax = plt.subplots(dpi=200)

    # Scatter plot of [Fe/H] true vs delta [Fe/H]
    ax.plot(feh_true, d_feh, marker='o', lw=0, color='black', markersize=4)
    # Linear fit to the data
    p = np.polyfit(feh_true, d_feh, deg=1)
    xx = np.linspace(-0.5, 0.5, 200)
    yy = np.polyval(p, xx)
    ax.plot(xx, yy, color='crimson', label='Linear fit\n$p_1$={:0.4f}\n$p_2$={:0.4f}'.format(p[0], p[1]))

    # Customize legend with specified font size
    ax.legend(loc='upper right', fontsize=legend_fontsize)

    # Set axis labels with custom font sizes
    ax.set_xlabel('$\mathrm{[Fe/H]}_{\mathrm{Archive}}$', fontsize=xlabel_fontsize)
    ax.set_ylabel('$\mathrm{[Fe/H]}_{\mathrm{Recovered}}$ - $\mathrm{[Fe/H]}_{\mathrm{Archive}}$', fontsize=ylabel_fontsize)

    # Customize tick labels and minor ticks with specified font sizes
    ax.tick_params(axis='both', which='major', labelsize=ticklabel_fontsize)
    ax.minorticks_on()


def plot_crossvalidation_results_main(order, df_crossval, savefolder=None,
                                      tick_fontsize=8, ticklabel_fontsize=9,
                                      xlabel_fontsize=10, ylabel_fontsize=10,
                                      legend_fontsize=7):
    """
    Main crossvalidation results plot with seaborn styling and proper layout adjustments.
    """
    PW = 10.
    PH = 8  # Adjust height to prevent overlapping
    fig = plt.figure(figsize=(PW, PH), dpi=300)
    gs = GridSpec(10, 10)
    gs.update(wspace=0.25, hspace=0.)

    ax = plt.subplot(gs[0:4, 0:3])
    bx = plt.subplot(gs[0:4, 3:6])
    cx = plt.subplot(gs[0:4, 7:])
    dx = plt.subplot(gs[5:, :])

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

    # Setting axis limits and labels with customized font sizes
    ax.set_xlim(-0.6, 0.6)
    ax.set_xlabel('[Fe/H]', labelpad=5, fontsize=xlabel_fontsize)
    ax.set_ylabel('$T_{\mathrm{eff}}$ [K]', labelpad=5, fontsize=ylabel_fontsize)
    bx.set_xlabel('$\log (g)$', labelpad=5, fontsize=xlabel_fontsize)
    ax.tick_params(axis='both', which='major', labelsize=ticklabel_fontsize)
    bx.tick_params(axis='both', which='major', labelsize=ticklabel_fontsize)
    ax.minorticks_on()

    # Remove y-ticks for the middle plot to avoid clutter
    bx.set_yticklabels([])
    bx.set_ylim(*ax.get_ylim())

    fig.subplots_adjust(wspace=0.15)

    # idx_cool = df_crossval['teff_true'].values < 4500
    # print(idx_cool)
    # Call the delta FeH plot function
    plot_crossval_feh_delta_feh(df_crossval.feh_true.values, df_crossval.d_feh.values, tick_fontsize=tick_fontsize,
                                ticklabel_fontsize=ticklabel_fontsize, xlabel_fontsize=xlabel_fontsize,
                                ylabel_fontsize=ylabel_fontsize, legend_fontsize=legend_fontsize, ax=cx)

    # Load data from pickle file
    with open(
            '/Users/tehan/Documents/SURFSUP/NEID_Spectra/NEIDSM/TIC 437229644_102/TIC 437229644_compositecomparison.pkl',
            'rb') as f:
        plot_data = pickle.load(f)

    # Extract data
    w = plot_data['wavelength']
    target_spectrum = plot_data['target_spectrum']
    composite_spectrum = plot_data['composite_spectrum']
    ref_spectra = plot_data['ref_spectra']
    labels = plot_data['labels']
    title = plot_data['title']
    scaleres = plot_data['scaleres']

    # Create plot for spectra
    dx.plot(w, target_spectrum, color='black', label='Target', lw=1)
    dx.plot(w, composite_spectrum, color='crimson', label='Composite', alpha=0.5, lw=1.5)

    for i, ref_spectrum in enumerate(ref_spectra):
        dx.plot(w, ref_spectrum + (5.0 - i), lw=1, color='black')
        dx.text(w[0], (6.2 - i), labels[i].replace('Teff', r'$T_{\text{eff}}$'), color='black', fontsize=tick_fontsize)

    dx.text(w[0], 1.2, 'Target Spectrum (Black), Composite Spectrum (Red)', fontsize=tick_fontsize)
    dx.text(w[0], 0.2, f'Residual = Target - Composite (RMS: {np.std(target_spectrum - composite_spectrum) * 1e3:.1f} ppt)', fontsize=tick_fontsize)
    print(np.sqrt(np.mean((target_spectrum - composite_spectrum)**2)))
    dx.set_title(title, fontsize=xlabel_fontsize)
    dx.plot(w, (target_spectrum - composite_spectrum) * scaleres, color='black', lw=1)
    dx.set_xlabel(r'Wavelength [$\AA$]', fontsize=xlabel_fontsize, labelpad=2)
    dx.set_ylabel('Normalize Flux + offset', fontsize=ylabel_fontsize, labelpad=2)
    dx.set_yticklabels([])
    dx.set_ylim(-1, 7)
    # dx.set_xlim(8655, 8665)

    # Customize tick font sizes
    dx.tick_params(axis='both', which='major', labelsize=ticklabel_fontsize)

    # Customize legend font size
    dx.legend(fontsize=legend_fontsize, loc='upper right')

    # Save or display the figure
    if savefolder:
        fig.savefig(f'{savefolder}/crossvalidation_o{order}_main.pdf', bbox_inches='tight')
    else:
        plt.show()


if __name__ == '__main__':
    order = 102
    df_crossval = pd.read_csv(
        f'/Users/tehan/PycharmProjects/neidspecmatch/library/20250226_specmatch_nir/crossval/o{order}_crossval/crossvalidation_results_o{order}.csv')
    print(np.where(np.abs(df_crossval['d_feh']) > 0.3))
    print(df_crossval.iloc[31])
    # Running the function with the provided data
    plot_crossvalidation_results_main(order=order, df_crossval=df_crossval,
                                      savefolder=f'/Users/tehan/Documents/SURFSUP/NEID_Spectra/NEIDSM/'
                                      )
