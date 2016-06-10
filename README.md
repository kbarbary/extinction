extinction
==========

*Fast interstellar dust extinction laws in Python*

[![Build Status](http://img.shields.io/travis/kbarbary/extinction.svg?style=flat-square&label=build)](https://travis-ci.org/kbarbary/extinction)

This module contains Cython-optimized implementations of a few
empirical dust exitinction laws found in the literature:

- [Cardelli, Clayton and Mathis (1989)](http://adsabs.harvard.edu/abs/1989ApJ...345..245C)
- [O'Donnell (1994)](http://adsabs.harvard.edu/abs/1994ApJ...422..158O)
- [Fitzpatrick (1999)](http://adsabs.harvard.edu/abs/1999PASP..111...63F)

## Install

```
pip install extinction
```

numpy and scipy are dependencies.


## API

Each function accepts a 1-d numpy array (not a list) of wavelengths and
some parameters, and returns the extinction in magnitudes at each wavelength.
The parameters are typically A_V, the total V band extinction in magnitudes
and R_V. R_V affects the shape of the returned function, while A_V only
affects its amplitude.

```python
import numpy as np
import extinction

wave = np.array([2000., 3000., 4000.])  # wavelength in Angstroms
a_v = 1.0
r_v = 3.1

# Cardelli, Clayton & Mathis (1989)
extinction.ccm89(wave, a_v, r_v)

# O'Donnell (1994)
extinction.odonnell94(wave, a_v, r_v)

# Fitzpatrick (1999) for R_V = 3.1
extinction.fitzpatrick99(wave, a_v)


# Fitzpatrick (1999) is based on splines, which have to be
# constructed for a given R_V. To use Fitzpatrick with R_V values other
# than 3.1:
f = extinction.Fitzpatrick99(r_v)  # construct once
f(wave, a_v)  # call multiple times

```

The above functions accept wavelength in Angstroms.
that accept wavenumber in inverse microns:

```python
x = 1e4 / np.array([2000., 3000., 4000.])  # wavenumber in inverse microns

extinction.ccm89_invum(x, a_v, r_v)
extinction.od94_invum(x, a_v, r_v)

f = extinction.F99Extinction(r_v)
f(x, a_v, unit='invum')
```

## Comparison of functions

![comparison plot](extinction.png)


## License

is MIT.