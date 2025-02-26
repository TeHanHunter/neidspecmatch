import matplotlib.pyplot as plt
import glob2
from tqdm import trange
import os
from astropy.io import ascii
import numpy as np
from astropy.io import fits


def plot_orders():
    filenames = glob2.glob('/home/tehan/Documents/hpf/spectra/*.fits')
    star_num = 0
    vsinis = np.zeros((11, len(filenames)))
    teffs = np.zeros((11, len(filenames)))
    logg = np.zeros((11, len(filenames)))
    feh = np.zeros((11, len(filenames)))
    star_id = []
    for i in trange(len(filenames)):
        target = 'TIC_' + os.path.basename(filenames[i]).split('_')[2]
        try:
            data = ascii.read(f'/home/tehan/Documents/hpf/{target}/{target}_overview.csv', format='csv')
            star_num += 1
            star_id.append(target)
        except FileNotFoundError:
            continue
        teffs[:, i] = data['teff']
        vsinis[:, i] = data['vsini']
        logg[:, i] = data['logg']
        feh[:, i] = data['feh']
        # median_teff = np.median(data['teff'])
        # median_vsini = np.median(data['vsini'])
        # colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'k']
        # for j in range(len(data)):
        #     order = data['filenames'][j].split('/')[-2].split('_')[-1]
        #     if star_num == 1:
        #         plt.scatter(star_num, np.log10(data['vsini'][j] / median_vsini), color=colors[j], marker='x',
        #                     linewidths=3, s=100, label=order)
        #     else:
        #         plt.scatter(star_num, np.log10(data['vsini'][j] / median_vsini), color=colors[j], marker='x',
        #                     linewidths=3, s=100)
    idx = np.argwhere(np.all(vsinis[..., :] == 0, axis=0))
    vsinis = np.delete(vsinis, idx, axis=1)
    teffs = np.delete(teffs, idx, axis=1)
    logg = np.delete(logg, idx, axis=1)
    feh = np.delete(feh, idx, axis=1)

    sort_id = np.argsort(np.median(teffs, axis=0))
    star_id = np.array(star_id)[sort_id]
    vsinis = vsinis[:, sort_id]
    teffs = teffs[:, sort_id]
    logg = logg[:, sort_id]
    feh = feh[:, sort_id]
    colors = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9', 'k']
    orders = ['14', '15', '16', '17', '18', '19', '26', '3', '4', '5', '6']

    fig = plt.figure(constrained_layout=False, figsize=(8, 10))
    gs = fig.add_gridspec(4, 1)
    gs.update(wspace=0.2, hspace=0.05)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[2, 0])
    ax4 = fig.add_subplot(gs[3, 0])

    for j in range(11):
        ax1.scatter(np.arange(len(vsinis[0])), np.log10(vsinis[j] / np.median(vsinis, axis=0)), color=colors[j],
                    marker='x',
                    linewidths=3, s=100, label=orders[j])
        ax1.plot(np.arange(len(vsinis[0])), np.log10(vsinis[j] / np.median(vsinis, axis=0)), c=colors[j], marker='')

        ax2.scatter(np.arange(len(vsinis[0])), (teffs[j] / np.median(teffs, axis=0)), color=colors[j], marker='x',
                    linewidths=3, s=100, label=orders[j])
        ax2.plot(np.arange(len(vsinis[0])), (teffs[j] / np.median(teffs, axis=0)), c=colors[j], marker='')

        ax3.scatter(np.arange(len(vsinis[0])), (logg[j] / np.median(logg, axis=0)), color=colors[j], marker='x',
                    linewidths=3, s=100, label=orders[j])
        ax3.plot(np.arange(len(vsinis[0])), (logg[j] / np.median(logg, axis=0)), c=colors[j], marker='')

        ax4.scatter(np.arange(len(vsinis[0])), (feh[j] / np.median(feh, axis=0)), color=colors[j], marker='x',
                    linewidths=3, s=100, label=orders[j])
        ax4.plot(np.arange(len(vsinis[0])), (feh[j] / np.median(feh, axis=0)), c=colors[j], marker='')

    # ax2.errorbar(np.arange(len(vsinis[0])), np.median(vsinis, axis=0), np.std(vsinis, axis=0), capsize=3, ls='')
    # ax3.errorbar(np.arange(len(vsinis[0])), np.median(teffs, axis=0), np.std(teffs, axis=0), capsize=3, ls='')
    for k in range(18):
        ax4.text(np.arange(len(vsinis[0]))[k] - 0.3, -55, star_id[k], verticalalignment='top', rotation=-90,
                 fontweight='semibold')
    # ax3.set_xlabel('Star No.')
    ax1.set_ylabel('log(vsini)')
    ax2.set_ylabel('teff')
    ax3.set_ylabel('logg')
    ax4.set_ylabel('feh')
    ax1.set_xticklabels([])
    ax2.set_xticklabels([])
    ax3.set_xticklabels([])
    ax4.set_xticklabels([])
    ax1.tick_params(axis='x', bottom=False)
    ax2.tick_params(axis='x', bottom=False)
    ax3.tick_params(axis='x', bottom=False)
    ax4.tick_params(axis='x', bottom=False)

    plt.gcf().subplots_adjust(bottom=0.15, top=0.95)
    ax4.legend(ncol=3)
    plt.savefig('/home/tehan/Documents/hpf/hpf_orders.pdf', dpi=300)
    plt.show()


def plot_lc_each_sector(local_directory=None):
    files = glob2.glob(f'{local_directory}*.fits')
    os.makedirs(f'{local_directory}plots/', exist_ok=True)
    for i in range(len(files)):
        with fits.open(files[i], mode='denywrite') as hdul:
            plt.figure(constrained_layout=False, figsize=(8, 4))
            plt.plot(hdul[1].data['time'], hdul[1].data['cal_psf_flux'], '.', label='cal_psf')
            plt.plot(hdul[1].data['time'], hdul[1].data['cal_aper_flux'], '.', label='cal_aper')
            plt.title(f'TIC_{hdul[0].header["TICID"]}')
            plt.legend()
            plt.savefig(f'{local_directory}plots/TIC_{hdul[0].header["TICID"]}.png', dpi=300)

def plot_lc(local_directory=None, catalog=None):
    targets = ascii.read(catalog, format='csv')
    os.makedirs(f'{local_directory}plots/', exist_ok=True)
    for i in trange(len(targets)):
        files = glob2.glob(f'{local_directory}*{targets["GaiaDR3"][i]}*.fits')
        if len(files) == 0:
            continue
        else:
            fig = plt.figure(figsize=(13, 5))
            for j in range(len(files)):
                with fits.open(files[j], mode='denywrite') as hdul:
                    q = list(hdul[1].data['TESS_flags'] == 0) and list(hdul[1].data['TGLC_flags'] == 0)
                    period = targets['period'][np.where(targets['TIC'] == int(hdul[0].header['TICID']))]
                    plt.plot(hdul[1].data['time'] % period, hdul[1].data['cal_aper_flux'], '.', c='silver', ms=1)
                    plt.plot(hdul[1].data['time'][q] % period, hdul[1].data['cal_aper_flux'][q], '.', c="C0", ms=1)
                    title = f'TIC_{hdul[0].header["TICID"]} with {len(files)} sector(s) of data'
            if targets['TIC'][i] == 236785891:
                plt.ylim(0.6, 1.4)
            plt.title(title)
            plt.xlabel('Phase (days)')
            plt.ylabel('Normalized flux')
            plt.savefig(f'{local_directory}/plots/{title}.png',
                        dpi=300)
            plt.close(fig)


import numpy as np
import matplotlib.pyplot as plt


def plot_vsini_limit(file='/Users/tehan/PycharmProjects/neidspecmatch/all_vsinis_he.npy'):
    all_vsinis = np.load(file, allow_pickle=True)
    injected_vsini = np.linspace(10, 0.01, 20)[:len(all_vsinis)]

    fig, axs = plt.subplots(2, 1, figsize=(6, 8), sharex=True, gridspec_kw={'height_ratios': [2, 1]})

    # Main plot: Injected vs. Recovered vsini
    axs[0].plot(injected_vsini, all_vsinis[:, 0], marker='.', label='i55')
    axs[0].plot(injected_vsini, all_vsinis[:, 1], marker='.', label='i101')
    axs[0].plot(injected_vsini, all_vsinis[:, 2], marker='.', label='i102')
    axs[0].plot(injected_vsini, all_vsinis[:, 3], marker='.', label='i103')
    axs[0].plot([0, 10], [0, 10], 'k')
    axs[0].set_ylabel('Recovered vsini')
    axs[0].legend()

    # Residual plot: Difference between recovered and injected
    for i in range(4):
        axs[1].plot(injected_vsini, all_vsinis[:, i] - injected_vsini, marker='.', label=f'Residual {i}')

    axs[1].axhline(0, color='k', linestyle='--')
    axs[1].set_xlabel('Injected vsini')
    axs[1].set_ylabel('Residual')

    plt.savefig('HE_vsini_recovery.png')
    plt.show()


if __name__ == '__main__':
    # plot_lc(local_directory = '/home/tehan/Documents/GEMS/lc/',
    #         catalog = '/home/tehan/Documents/GEMS/GEMS.csv')
    plot_vsini_limit()