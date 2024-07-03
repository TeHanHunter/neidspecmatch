import os
import numpy as np
from astropy.io import ascii, fits
from astropy import units as u
from astropy import constants
from astroquery.simbad import Simbad
from astropy.table import Table, Column
from astropy.coordinates import SkyCoord
from tqdm import trange
DIRNAME = os.path.dirname(os.path.dirname(__file__))

def load_Mann(file='asu.fit'):
    Mann = fits.open(f'{DIRNAME}/library/{file}')
    return Mann[1].data


def load_Yee(file='apjaa58eat6_mrt.txt'):
    Yee = ascii.read(f'{DIRNAME}/library/{file}')
    return Yee

#from Caleb
def calclogg(mass, masserr, radius, raderr):
    mass = np.array(mass, copy=True) * u.Msun
    masserr = np.array(masserr, copy=True) * u.Msun
    radius = np.array(radius, copy=True) * u.Rsun
    raderr = np.array(raderr, copy=True) * u.Rsun
    val = np.log10((constants.G * mass * radius ** (-2)).cgs.value)
    valerr = np.sqrt(
        ((mass ** (-2.) * masserr ** 2. + 4. * raderr ** 2. * radius ** (-2.)) * np.log(10.) ** (-2.)).value)
    return val, valerr

def combine_Mann_Yee():
    Mann = load_Mann()
    Yee = load_Yee()
    print(Mann[-20])
    combined = Yee
    combined.add_column(Column(name='Source', data=['Yee'] * len(Yee)))
    for i in range(len(Mann)):
        logg, e_logg = calclogg(Mann[i][51], Mann[i][52], Mann[i][49], Mann[i][50])
        combined.add_row({'Name': f'{Mann[i][3]}', 'Teff': f'{Mann[i][47]}', 'e_Teff': f'{Mann[i][48]}',
                                        'R*': f'{Mann[i][49]:.3f}', 'e_R*': f'{Mann[i][50]:.3f}', 'log(g)': f'{logg:.3f}',
                                        'e_log(g)': f'{e_logg:.3f}', '[Fe/H]': f'{Mann[i][14]:.3f}', 'e_[Fe/H]': f'{Mann[i][15]:.3f}',
                                        'M*': f'{Mann[i][51]:.3f}', 'e_M*': f'{Mann[i][52]:.3f}',
                                        'logA': f'{Mann[i][53]:.3f}', 'e_logA': f'{Mann[i][54]:.3f}', 'plx': 0.,
                                        'Vmag': f'{Mann[i][19]:.3f}', 'Notes': 'N/A', 'Source': 'Mann'})
    ra = []
    dec = []
    for j in trange(len(combined)):
        try:
            result_table = Simbad.query_object(combined['Name'][j])
            ra_deg = result_table['RA']
            dec_deg = result_table['DEC']
            coord = SkyCoord(ra=ra_deg[0], dec=dec_deg[0], unit=(u.hourangle, u.deg))
            ra.append(coord.ra.deg)
            dec.append(coord.dec.deg)
        except TypeError:
            print(j, combined[j]['Name'])
            ra.append(-100.)
            dec.append(-100.)
    combined.add_column(Column(name='RA', data=ra))
    combined.add_column(Column(name='Dec', data=dec))
    return combined

def find_star(name, combined=None):
    custom_simbad = Simbad()
    custom_simbad.add_votable_fields('rvz_radvel')
    result_table = custom_simbad.query_object(name)
    radial_velocity = result_table['RVZ_RADVEL'][0]
    print(f"Radial Velocity of {name}: {radial_velocity} km/s")
    ra_deg = result_table['RA']
    dec_deg = result_table['DEC']
    coord = SkyCoord(ra=ra_deg[0], dec=dec_deg[0], unit=(u.hourangle, u.deg))
    idx = np.argmin((combined['RA']-coord.ra.deg)**2 + (combined['Dec']-coord.dec.deg)**2)
    print(np.sqrt((combined['RA']-coord.ra.deg)**2 + (combined['Dec']-coord.dec.deg)**2)[idx])
    print(combined[idx])

if __name__ == '__main__':
    name = 'HD 141004'
    # combined = combine_Mann_Yee()
    # combined.write(f'{DIRNAME}/library/combined.csv', format='csv',)
    combined = Table.read(f'{DIRNAME}/library/combined.csv', format='csv',)
    find_star(name, combined=combined)
    # All simbad IDs
    simbadres = Simbad.query_objectids(name)
    # simbadres.pprint_all()
    all_IDS = ''
    for i in range(len(simbadres)):
        all_IDS += str(simbadres[i][0])
        all_IDS += '|'
    print(all_IDS)
