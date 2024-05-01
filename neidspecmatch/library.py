import os
import sigfig
import numpy as np
from astropy.io import ascii, fits
from astropy import units as u
from astropy import constants
from astroquery.simbad import Simbad

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

def find_star(name='HD 100623'):
    Mann = load_Mann()
    Yee = load_Yee()
    Mann_names = []
    for i in range(len(Mann)):
        Mann_names.append(Mann[i][3])
    Mann_names = np.array(Mann_names)
    Yee_names = np.array(Yee['Name'])
    # print(Mann_names)
    # print(Yee_names[0])
    Mann_row = np.where(Mann_names == name)
    Yee_row = np.where(Yee_names == name)
    return Mann[Mann_row], Yee[Yee_row]

if __name__ == '__main__':
    name = 'HIP 23512'
    Mann_star, Yee_star = find_star(name=name)
    print(Mann_star)
    print(calclogg(0.44,0.044,0.443,0.017))
    Yee_star.pprint_all()

    # All simbad IDs
    simbadres = Simbad.query_objectids(name)
    # simbadres.pprint_all()
    all_IDS = ''
    for i in range(len(simbadres)):
        all_IDS += str(simbadres[i][0])
        all_IDS += '|'
    print(all_IDS)
