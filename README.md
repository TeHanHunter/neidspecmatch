# NEIDspecmatch (in progress)
NEIDSpecMatch: Spectral matching of NEID data

# Dependencies 

- pyde, either (`pip install pyde`) or install from here: https://github.com/hpparvi/PyDE This package needs numba (try `conda install numba` if problems).
- emcee (`pip install emcee`)
- astroquery (`pip install astroquery`)
- crosscorr (`git clone https://github.com/TeHanHunter/crosscorr.git`) `pip3 install .` NEED fortran installation. For Mac: brew install gcc (GNU fortran). For Ubuntu: sudo apt install gfortran
- NEIDspec (`git clone https://github.com/TeHanHunter/neidspec.git`) `pip3 install .`
- lmfit (`pip install lmfit`)
- barycorrpy (`pip install barycorrpy`)

# Installation
create a new conda env with
conda create -n neidspecmatch python==3.10
conda activate neidspecmatch
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
git clone https://github.com/TeHanHunter/neidspecmatch.git
cd neidspecmatch
pip3 install .
```
