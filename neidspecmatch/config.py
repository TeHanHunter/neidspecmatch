import os
import glob

# HPF wavelength bounds for different orders in Angstrom
BOUNDS = {'55': [5160., 5230.],
          '101': [8435., 8585.],
          '102': [8565., 8700.],
          '103': [8680., 8835.],}

# Directory name of package
DIRNAME = os.path.dirname(os.path.dirname(__file__))
# DIRNAME = '/home/sejones/neidspecmatch'
print('DIRNAME: {}'.format(DIRNAME))

# Default library path
PATH_LIBRARIES = os.path.join(DIRNAME, "library")
PATH_LIBRARY = os.path.join(PATH_LIBRARIES, "20240301_specmatch_nir")
PATH_LIBRARY_DB = os.path.join(PATH_LIBRARY, "20240301_46stars.csv")
PATH_LIBRARY_FITS = os.path.join(PATH_LIBRARY, "FITS")
PATH_LIBRARY_CROSSVAL = os.path.join(PATH_LIBRARY, "crossval")
PATH_LIBRARY_ZIPNAME = os.path.join(PATH_LIBRARY, '20210406_specmatch_nir_library.zip')
URL_LIBRARY = 'https://www.dropbox.com/s/rtees0v6yt0t9eb/20210811_specmatch_nir.zip?dl=1'
LIBRARY_FITSFILES = sorted(glob.glob(PATH_LIBRARY_FITS + '/*.fits'))
print(PATH_LIBRARY_FITS)
