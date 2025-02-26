import numpy as np
import matplotlib.pylab as plt
from PyAstronomy import pyasl

# Create data with a Gaussian absoprtion line
wvl = np.arange(4999., 5011., 0.04)
flux = np.zeros(len(wvl))

# The Gaussian
A = -0.05
s = 0.1
mu = 5004.1635788
flux += A/np.sqrt(2.*np.pi*s**2) * \
    np.exp(-(wvl-mu)**2/(2.*s**2))

# Apply the fast algorithm and ...
bfast = pyasl.fastRotBroad(wvl, flux, 0.81, 20)
# ... the slower one
bslow = pyasl.rotBroad(wvl, flux, 0.81, 20)

if __name__ == '__main__':
    plt.xlabel("Wvl [A]")
    plt.ylabel("Flux [au]")
    plt.title("Initial spectrum (black), fast (blue), slow (red, shifted)")
    plt.plot(wvl, flux, 'k.-')
    plt.plot(wvl, bfast, 'b.-')
    plt.plot(wvl, bslow, 'r.-')
    plt.show()
