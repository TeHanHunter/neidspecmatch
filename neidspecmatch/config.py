import os
import glob

# HPF wavelength bounds for different orders in Angstrom
BOUNDS = {'10': [3734, 3786],
          '11': [3757, 3810],
          '12': [3781, 3833],
          '13': [3804, 3858],
          '14': [3828, 3882],
          '15': [3852, 3906],
          '16': [3877, 3931],
          '17': [3901, 3957],
          '18': [3927, 3982],
          '19': [3952, 4008],
          '20': [3978, 4034],
          '21': [4004, 4061],
          '22': [4030, 4088],
          '23': [4057, 4115],
          '24': [4084, 4143],
          '25': [4112, 4171],
          '26': [4140, 4200],
          '27': [4168, 4228],
          '28': [4197, 4258],
          '29': [4226, 4287],
          '30': [4255, 4317],
          '31': [4285, 4348],
          '32': [4316, 4379],
          '33': [4346, 4410],
          '34': [4377, 4442],
          '35': [4409, 4474],
          '36': [4441, 4507],
          '37': [4474, 4540],
          '38': [4507, 4574],
          '39': [4540, 4608],
          '40': [4575, 4643],
          '41': [4609, 4678],
          '42': [4644, 4714],
          '43': [4680, 4750],
          '44': [4716, 4787],
          '45': [4753, 4824],
          '46': [4790, 4862],
          '47': [4828, 4901],
          '48': [4867, 4940],
          '49': [4906, 4980],
          '50': [4946, 5021],
          '51': [4986, 5062],
          '52': [5027, 5104],
          '53': [5069, 5147],
          '54': [5112, 5190],
          '55': [5155, 5234],
          '56': [5199, 5279],
          '57': [5244, 5324],
          '58': [5289, 5371],
          '59': [5335, 5418],
          '60': [5383, 5466],
          '61': [5431, 5515],
          '62': [5479, 5565],
          '63': [5529, 5615],
          '64': [5580, 5667],
          '65': [5631, 5720],
          '66': [5684, 5773],
          '67': [5737, 5828],
          '68': [5792, 5883],
          '69': [5848, 5940],
          '70': [5904, 5998],
          '71': [5962, 6057],
          '72': [6021, 6117],
          '73': [6081, 6178],
          '74': [6142, 6240],
          '75': [6205, 6304],
          '76': [6269, 6369],
          '77': [6334, 6436],
          '78': [6401, 6504],
          '79': [6469, 6573],
          '80': [6538, 6644],
          '81': [6609, 6716],
          '82': [6682, 6790],
          '83': [6756, 6865],
          '84': [6831, 6943],
          '85': [6909, 7022],
          '86': [6988, 7103],
          '87': [7069, 7185],
          '88': [7153, 7270],
          '89': [7238, 7357],
          '90': [7325, 7445],
          '91': [7414, 7536],
          '92': [7505, 7629],
          '93': [7599, 7725],
          '94': [7695, 7823],
          '95': [7794, 7923],
          '96': [7895, 8026],
          '97': [7998, 8132],
          '98': [8105, 8241],
          '99': [8214, 8352],
          '100': [8327, 8467],
          '101': [8442, 8584],
          '102': [8561, 8705],
          '103': [8683, 8830]}

# Directory name of package
DIRNAME = os.path.dirname(os.path.dirname(__file__))
# DIRNAME = '/home/sejones/neidspecmatch'
print('DIRNAME: {}'.format(DIRNAME))

# Default library path
PATH_LIBRARIES = os.path.join(DIRNAME, "library")
PATH_LIBRARY = os.path.join(PATH_LIBRARIES, "20240822_specmatch_nir")
PATH_LIBRARY_DB = os.path.join(PATH_LIBRARY, "20240822_89stars.csv")
PATH_LIBRARY_FITS = os.path.join(PATH_LIBRARY, "FITS")
PATH_LIBRARY_CROSSVAL = os.path.join(PATH_LIBRARY, "crossval")
PATH_LIBRARY_ZIPNAME = os.path.join(PATH_LIBRARY, '20210406_specmatch_nir_library.zip')
URL_LIBRARY = 'https://www.dropbox.com/s/rtees0v6yt0t9eb/20210811_specmatch_nir.zip?dl=1'
LIBRARY_FITSFILES = sorted(glob.glob(PATH_LIBRARY_FITS + '/*.fits'))
print(PATH_LIBRARY_FITS)
