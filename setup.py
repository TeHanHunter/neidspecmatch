from setuptools import setup, find_packages

def readme():
    with open('README.md') as f:
        return f.read()

setup(
    name='neidspecmatch',
    version='0.1.3',
    description='Matching NEID Spectra',
    long_description=readme(),
    long_description_content_type='text/markdown',
    url='https://github.com/gummiks/hpfspecmatch/',
    author='Te Han, Gudmundur Stefansson (author of HPFSpecMatch)',
    author_email='tehanhunter@gmail.com',
    install_requires=[
        'numpy == 1.24.3',
        'numba',
        'pytransit',
        'numpy',
        'barycorrpy',
        'emcee',
        'lmfit',
        'neidspec',
        'crosscorr',
        'astroquery==0.4.7',
        'glob2',
        'celerite',
        'radvel',
        'wget',
        'siphash24',
    ],
    packages=find_packages(),
    license_files='MIT',
    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Topic :: Scientific/Engineering :: Astronomy',
    ],
    keywords='NEID Spectra Astronomy',
    include_package_data=True,
    zip_safe=False,
)