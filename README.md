# NEIDspecmatch (in progress)
NEIDSpecMatch: Spectral matching of NEID data. Modified from [HPFSpecMatch](https://gummiks.github.io/hpfspecmatch/). 

# Dependencies 

- pyde, either (`pip install pyde`) or install from [here](https://github.com/hpparvi/PyDE). This package needs numba (try `conda install numba` if problems).
- emcee (`pip install emcee`)
- astroquery (`pip install astroquery`)
- crosscorr (`git clone https://github.com/TeHanHunter/crosscorr.git`) `pip3 install .` NEED fortran installation. For Mac: brew install gcc (GNU fortran). For Ubuntu: sudo apt install gfortran
- NEIDspec (`git clone https://github.com/TeHanHunter/neidspec.git`) `pip3 install .`
- lmfit (`pip install lmfit`)
- barycorrpy (`pip install barycorrpy`)

  Known Issue: The latest version of barycorrpy deprecated some syntax used in the NEIDSpecMatch. Please use earlier versions (0.4.4 tested to work) while we update the syntax. 

# Installation
create a new conda env with
`conda create -n neidspecmatch python==3.10`
`conda activate neidspecmatch`
```
conda install numba
git clone https://github.com/hpparvi/PyDE.git
cd PyDE
pip3 install .
cd ..
pip3 install emcee
pip3 install astroquery
git clone https://github.com/TeHanHunter/crosscorr.git
cd crosscorr
brew install gcc
pip3 install .
cd ..
git clone https://github.com/TeHanHunter/neidspec.git
cd neidspec
pip3 install .
cd ..
pip3 install lmfit
pip3 install barycorrpy
pip3 install celerite
git clone https://github.com/TeHanHunter/neidspecmatch.git
cd neidspecmatch
pip3 install .
```
# Library download
The current library include 89 library stars. You may download them [here](https://drive.google.com/drive/folders/1_ZcQavq5boQt5f7RjelEUF6sqZKibSbT?usp=sharing).

# Cross-validation and result uncertainty 
The cross-validation is a necessary step to estimate the uncertainty of the recovered stellar parameters. To run the cross-validation on the current library, one can run `run_crossval.py` for the desired order. The output include a file in the format of `crossvalidation_results_o102.csv`. One can take the standard deviation of the column d_teff, d_feh, and	d_logg on certain rows (depending on the star type you want to estimate the uncertainty). 

# Fit a spectrum
Once the library is in place, one can run NEIDSpecMatch using `run_neidspecmatch.py`. 

# Reference
TODO
