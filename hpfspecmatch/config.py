import os
import glob
# HPF wavelength bounds for different orders in Angstrom
BOUNDS = {'0': [8085., 8175.], 
          '1': [8195., 9285.],
          '2': [8310., 8400.],
          '3': [8425., 8515.],
          ### original
          '4': [8540.,8640.],
          '5': [8670.,8750.],
          '6': [8790.,8885.],
          ###
          '7': [8920., 9015.],
          '8': [9050., 9150.],
          '9': [9190., 9290.],
          '10': [9330., 9435.],
          '11': [9475., 9580.],
          '12': [9625., 9735.],
          '13': [9780., 9890.],
          ### original
          '14': [9940.,10055.],
          '15': [10105.,10220.],
          '16': [10280.,10395.],
          '17': [10460.,10570.],
          ###
          '18': [10640., 10760.],
          '19': [10830., 10950.],
          '20': [11025., 11150.],
          '21': [11230., 11355.],
          '22': [11440., 11570.],
          '23': [11660., 11790.],
          '24': [11890., 12025.],
          '25': [12125., 12265.],
          '26': [12375., 12520.],
          '27': [12630., 12780.]}

# Directory name of package
DIRNAME = os.path.dirname(os.path.dirname(__file__))
#DIRNAME = '/home/sejones/hpfspecmatch'
print('DIRNAME: {}'.format(DIRNAME))

# Default library path
PATH_LIBRARIES = os.path.join(DIRNAME,"library")
#PATH_LIBRARY = os.path.join(PATH_LIBRARIES,"20210406_specmatch_nir_library")
PATH_LIBRARY = os.path.join(PATH_LIBRARIES,"20210811_specmatch_nir")
PATH_LIBRARY_DB = os.path.join(PATH_LIBRARY,"20210808_spec_mann_w_files_qual1_pass3.csv")
PATH_LIBRARY_FITS = os.path.join(PATH_LIBRARY,"FITS")
PATH_LIBRARY_CROSSVAL = os.path.join(PATH_LIBRARY,"crossval")
PATH_LIBRARY_ZIPNAME = os.path.join(PATH_LIBRARY,'20210406_specmatch_nir_library.zip')
#URL_LIBRARY = 'https://www.dropbox.com/s/69j00zrpov48qwx/20201008_specmatch_nir.zip?dl=1'
#URL_LIBRARY = 'https://www.dropbox.com/sh/8vo1hgn118irjfy/AADaHzEP9XMqTgKlkNQjVTWWa?dl=1'
URL_LIBRARY = 'https://www.dropbox.com/s/rtees0v6yt0t9eb/20210811_specmatch_nir.zip?dl=1'
LIBRARY_FITSFILES = sorted(glob.glob(PATH_LIBRARY_FITS+'/*.fits'))
print(PATH_LIBRARY_FITS)
