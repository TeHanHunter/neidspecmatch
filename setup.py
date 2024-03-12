from setuptools import setup, find_packages

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='neidspecmatch',
      version='0.1.1',
      description='Matching HPF Spectra',
      long_description=readme(),
      url='https://github.com/gummiks/hpfspecmatch/',
      author='Gudmundur Stefansson',
      author_email='gummiks@gmail.com',
      install_requires=['barycorrpy','emcee','lmfit','neidspec','crosscorr','pyde','astroquery','glob2'],
      packages=['neidspecmatch'],
      license='GPLv3',
      classifiers=['Topic :: Scientific/Engineering :: Astronomy'],
      keywords='HPF Spectra Astronomy',
      dependency_links=['http://github.com/user/repo/tarball/master#egg=package-1.0'],
      include_package_data=True,
      zip_safe=False
      )
