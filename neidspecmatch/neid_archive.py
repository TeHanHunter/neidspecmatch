import sys
import os
import numpy as np
from pyneid.neid import Neid
from astropy.time import Time
import astropy.units as u
from astropy.table import Table, Column
from astropy.coordinates import SkyCoord
from tqdm import trange

DIRNAME = os.path.dirname(os.path.dirname(__file__))


def get_table_singlestar(starname, ra, dec, fmt, download=False, directory=None):
    print(starname)

    param = dict()
    param['datalevel'] = 'l2'
    param['object'] = starname
    param['radius'] = 1

    """
    other param options include:
        piname (e.g., "Suvrath Mahadevan")
        program (e.g., "2021B-2015")
        datetime (e.g., str(Time.now()-30*u.day)[:10]+" 00:00:00/")
    """

    if fmt == 'csv':
        outpath = starname.replace(' ', '_') + '_L2_q.csv'
    elif fmt == 'ipac':
        outpath = starname.replace(' ', '_') + '_L2_q.tbl'

    # Neid.query_criteria(param, \
    #                     cookiepath=f'{directory}neidadmincookie.txt', \
    #                     format=fmt, \
    #                     outpath=directory + outpath)
    Neid.query_position('l2', f'circle {ra} {dec} 0.5', cookiepath=f'{directory}neidadmincookie.txt', format=fmt,
                        outpath=directory + outpath)
    if fmt == 'ipac' and download == True:
        Neid.download(outpath,
                      'l2', \
                      'ipac', \
                      '.', \
                      cookiepath=f'{directory}neidadmincookie.txt')


if __name__ == '__main__':
    directory = '/Users/tehan/Documents/NEID_archive/'
    Neid.login(userid='than', password='QWAqwab-6766', cookiepath=f'{directory}neidadmincookie.txt',
               debugfile=f'{directory}archive.debug')
    combined = Table.read(f'{DIRNAME}/library/combined.csv', format='csv', )
    have_spectra = []
    for i in trange(len(combined)):  # len(combined)
        star = combined['Name'][i]
        ra = combined['RA'][i]
        dec = combined['Dec'][i]
        try:
            get_table_singlestar(star, ra, dec, 'csv', download=False, directory=directory)
            t = Table.read(f'{directory}{star.replace(" ", "_")}_L2_q.csv', format='csv')
            if len(t) > 0:
                have_spectra.append(star)
        except:
            continue
    print(have_spectra)
    result = Table([have_spectra], names=['Name'])
    result.write(f'{directory}have_spectra.csv', format='csv', overwrite=True)

    """
    if 'download' in sys.argv:
        fmt = 'ipac'
        download = True
    else:
        download = False
        fmt = 'csv'
    Example
    
    python3 -m pyneid_example_te user pw HD\ 143761
    
    then
    
    python3 -m pyneid_example_te user pw HD\ 143761 download
    """
