# NEIDspecmatch (in progress)
NEIDSpecMatch: Spectral matching of NEID data

# Dependencies 

- pyde, either (`pip install pyde`) or install from here: https://github.com/hpparvi/PyDE This package needs numba (try `conda install numba` if problems).
- emcee (`pip install emcee`)
- astroquery (`pip install astroquery`)
- crosscorr (`git clone https://github.com/TeHanHunter/crosscorr.git`) `pip3 install .`
- NEIDspec (`git clone https://github.com/TeHanHunter/neidspec.git`) `pip3 install .`
- lmfit (`pip install lmfit`)
- barycorrpy (`pip install barycorrpy`)

# Installation
```
conda install numba
pip install pyde
pip install emcee
pip install astroquery
git clone https://github.com/TeHanHunter/crosscorr.git
cd crosscorr
pip3 install .
git clone https://github.com/TeHanHunter/neidspec.git
cd neidspec
pip3 install .
pip install lmfit
pip install barycorrpy
git clone https://github.com/TeHanHunter/neidspecmatch.git
cd neidspecmatch
pip3 install .
```
