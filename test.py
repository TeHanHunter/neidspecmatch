import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import shutil
import os

if __name__ == '__main__':
    df_lib = pd.read_csv('/Users/tehan/PycharmProjects/neidspecmatch/library/20250204_specmatch_nir/20250204_89stars.csv')
    print(len(set(df_lib['qprog'])))
    data = df_lib['[Fe/H]']

    mu, sigma = norm.fit(data)
    print(mu, sigma)
    plt.hist(data, bins=15, density=True, alpha=0.6, color='k')
    x = np.linspace(data.min(), data.max(), 100)
    plt.plot(x, norm.pdf(x, mu, sigma), 'r-', lw=2)
    plt.xlabel('[Fe/H]')
    plt.show()

    # Generate samples with boundary constraints
    samples = []
    while len(samples) < len(data):
        sample = norm.rvs(mu, sigma)
        if -0.5 <= sample <= 0.5:
            samples.append(sample)

    plt.hist(data - samples, bins=15, density=True, alpha=0.6, color='k')
    print(np.std(data - samples))
    plt.xlabel('[Fe/H]')
    plt.show()

    filtered_df = df_lib[df_lib['Teff'] > 4000]
    source_dir = '/Users/tehan/PycharmProjects/neidspecmatch/library/20250204_specmatch_nir/FITS'
    destination_dir = '/Users/tehan/PycharmProjects/neidspecmatch/library/20250204_specmatch_nir/FITS_solar'

    # Move files
    for basename in filtered_df['basenames']:
        source_path = os.path.join(source_dir, basename)
        destination_path = os.path.join(destination_dir, basename)

        # Check if the file exists in the source directory
        if os.path.exists(source_path):
            shutil.move(source_path, destination_path)
            print(f"Moved: {basename}")
        else:
            print(f"File not found: {basename}")

    print(len(data))
    print(samples)
    # Replace only the relevant rows in the 'Fe/H' column
    df_lib.iloc[:len(samples), df_lib.columns.get_loc('[Fe/H]')] = samples
    df_lib = df_lib[df_lib['Teff'] <= 4000]
    df_lib.to_csv('20250204_m_dwarf_random_met.csv', index=False)