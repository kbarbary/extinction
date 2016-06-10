#!python
#cython: boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True

"""Interstellar dust extinction functions."""

import numpy as np
cimport numpy as np
from scipy.interpolate import splmake, spleval

__version__ = "0.2.0"

__all__ = ['ccm89', 'odonnell94', 'Fitzpatrick99', 'fitzpatrick99', 'fm07',
           'calzetti00']

# ------------------------------------------------------------------------------
# Utility functions for converting wavelength units

cdef double aa_to_invum(double x):
    """Convert Angstroms to inverse microns"""
    return 1e4 / x

cdef double noop(double x):
    return x

ctypedef double (*scalar_func)(double)

# -----------------------------------------------------------------------------
# Cardelli, Clayton & Mathis (1989)

cdef inline void ccm89ab_ir_invum(double x, double *a, double *b):
    """ccm89 a, b parameters for 0.3 < x < 1.1 (infrared)"""
    cdef double y
    y = x**1.61
    a[0] = 0.574 * y
    b[0] = -0.527 * y


cdef inline void ccm89ab_opt_invum(double x, double *a, double *b):
    """ccm89 a, b parameters for 1.1 < x < 3.3 (optical)"""

    cdef double y

    y = x - 1.82
    a[0] = ((((((0.329990*y - 0.77530)*y + 0.01979)*y + 0.72085)*y - 0.02427)*y
             - 0.50447)*y + 0.17699)*y + 1.0
    b[0] = ((((((-2.09002*y + 5.30260)*y - 0.62251)*y - 5.38434)*y + 1.07233)*y
             + 2.28305)*y + 1.41338)*y


cdef inline void ccm89ab_uv_invum(double x, double *a, double *b):
    """ccm89 a, b parameters for 3.3 < x < 8.0 (ultraviolet)"""
    cdef double y, y2, y3

    y = x - 4.67
    a[0] = 1.752 - 0.316*x - (0.104 / (y*y + 0.341))
    y = x - 4.62
    b[0] = -3.090 + 1.825*x + (1.206 / (y*y + 0.263))
    if x > 5.9:
        y = x - 5.9
        y2 = y * y
        y3 = y2 * y
        a[0] += -0.04473*y2 - 0.009779*y3
        b[0] += 0.2130*y2 + 0.1207*y3


cdef inline void ccm89ab_fuv_invum(double x, double *a, double *b):
    """ccm89 a, b parameters for 8 < x < 11 (far-UV)"""
    cdef double y, y2, y3

    y = x - 8.
    y2 = y * y
    y3 = y2 * y
    a[0] = -0.070*y3 + 0.137*y2 - 0.628*y - 1.073
    b[0] = 0.374*y3 - 0.420*y2 + 4.257*y + 13.670


cdef inline void ccm89ab_invum(double x, double *a, double *b):
    """ccm89 a, b parameters for 0.3 < x < 11 in microns^-1"""
    if x < 1.1:
        ccm89ab_ir_invum(x, a, b)
    elif x < 3.3:
        ccm89ab_opt_invum(x, a, b)
    elif x < 8.:
        ccm89ab_uv_invum(x, a, b)
    else:
        ccm89ab_fuv_invum(x, a, b)


def ccm89(double[:] wave, double a_v, double r_v, unit='aa',
          np.ndarray out=None):
    """ccm89(wave, a_v, r_v, unit='aa', out=None)

    Cardelli, Clayton & Mathis (1989) extinction function.

    Parameters
    ----------
    wave : numpy.ndarray (1-d)
        Wavelengths or wavenumbers.
    a_v : float
        Scaling parameter, A_V: extinction in magnitudes at characteristic
        V band wavelength.
    r_v : float
        Ratio of total to selective extinction, A_V / E(B-V).
    unit : {'aa', 'invum'}, optional
        Unit of wave: 'aa' (Angstroms) or 'invum' (inverse microns).
    out : np.ndarray, optional
        If specified, store output values in this array.

    Returns
    -------
    Extinction in magnitudes at each input wavelength.
    """

    cdef:
        size_t i
        size_t n = wave.shape[0]
        double a = 0.0
        double b = 0.0

    if out is None:
        out = np.empty(n, dtype=np.float)
    else:
        assert out.shape == wave.shape
        assert out.dtype == np.float

    cdef scalar_func convert_wave
    if unit == 'aa':
        convert_wave = &aa_to_invum
    elif unit == 'invum':
        convert_wave = &noop
    else:
        raise ValueError("unrecognized unit")

    cdef double[:] out_view = out
    for i in range(n):
        ccm89ab_invum(convert_wave(wave[i]), &a, &b)
        out_view[i] = a_v * (a + b / r_v)

    return out


# -----------------------------------------------------------------------------
# O'Donnell (1994)

cdef inline void od94ab_opt_invum(double x, double *a, double *b):
    """od94 a, b parameters for 1.1 < x < 3.3 (optical)"""
    cdef double y

    y = x - 1.82
    a[0] = (((((((-0.505*y + 1.647)*y - 0.827)*y - 1.718)*y + 1.137)*y +
              0.701)*y - 0.609)*y + 0.104)*y + 1.0
    b[0] = (((((((3.347*y - 10.805)*y + 5.491)*y + 11.102)*y - 7.985)*y -
              3.989)*y + 2.908)*y + 1.952)*y


cdef inline void od94ab_invum(double x, double *a, double *b):
    """od94 a, b parameters for 0.3 < x < 11 in microns^-1"""
    if x < 1.1:
        ccm89ab_ir_invum(x, a, b)
    elif x < 3.3:
        od94ab_opt_invum(x, a, b)
    elif x < 8.:
        ccm89ab_uv_invum(x, a, b)
    else:
        ccm89ab_fuv_invum(x, a, b)


def odonnell94(double[:] wave, double a_v, double r_v, unit='aa',
         np.ndarray out=None):
    """odonnell94(wave, a_v, r_v, out=None, unit='aa', out=None)

    O'Donnell (1994) extinction function.

    Like Cardelli, Clayton, & Mathis (1989) [1]_ but using the O'Donnell
    (1994) [2]_ optical coefficients between 3030 A and 9091 A.

    Parameters
    ----------
    wave : numpy.ndarray (1-d)
        Wavelengths or wavenumbers.
    a_v : float
        Scaling parameter, A_V: extinction in magnitudes at characteristic
        V band wavelength.
    r_v : float
        Ratio of total to selective extinction, A_V / E(B-V).
    unit : {'aa', 'invum'}, optional
        Unit of wave: 'aa' (Angstroms) or 'invum' (inverse microns).
    out : np.ndarray, optional
        If specified, store output values in this array.

    Returns
    -------
    Extinction in magnitudes at each input wavelength.

    Notes
    -----
    This function matches the Goddard IDL astrolib routine CCM_UNRED.
    From the documentation for that routine:
    1. The CCM curve shows good agreement with the Savage & Mathis (1979)
       [3]_ ultraviolet curve shortward of 1400 A, but is probably
       preferable between 1200 and 1400 A.
    2. Curve is extrapolated between 912 and 1000 A as suggested by
       Longo et al. (1989) [4]_
    3. Valencic et al. (2004) [5]_ revise the ultraviolet CCM
       curve (3.3 -- 8.0 um^-1).    But since their revised curve does
       not connect smoothly with longer and shorter wavelengths, it is
       not included here.

    References
    ----------
    .. [1] Cardelli, J. A., Clayton, G. C., & Mathis, J. S. 1989, ApJ, 345, 245
    .. [2] O'Donnell, J. E. 1994, ApJ, 422, 158O 
    .. [3] Savage & Mathis 1979, ARA&A, 17, 73
    .. [4] Longo et al. 1989, ApJ, 339,474
    .. [5] Valencic et al. 2004, ApJ, 616, 912
    """

    cdef:
        size_t i
        size_t n = wave.shape[0]
        double a = 0.0
        double b = 0.0

    if out is None:
        out = np.empty(n, dtype=np.float)
    else:
        assert out.shape == wave.shape
        assert out.dtype == np.float

    cdef scalar_func convert_wave
    if unit == 'aa':
        convert_wave = &aa_to_invum
    elif unit == 'invum':
        convert_wave = &noop
    else:
        raise ValueError("unrecognized unit")

    cdef double[:] out_view = out
    for i in range(n):
        od94ab_invum(convert_wave(wave[i]), &a, &b)
        out_view[i] = a_v * (a + b / r_v)

    return out


# -----------------------------------------------------------------------------
# Fitzpatrick (1999)

DEF F99_X0 = 4.596
DEF F99_GAMMA = 0.99
DEF F99_C3 = 3.23
DEF F99_C4 = 0.41
DEF F99_C5 = 5.9
DEF F99_X02 = F99_X0 * F99_X0
DEF F99_GAMMA2 = F99_GAMMA * F99_GAMMA


cdef inline double f99_uv_invum(double x, double a_v, double r_v):
    """f99 extinction for x < 1e4/2700 microns^-1"""
    cdef double c1, c2, d, x2, y, y2, rv2, k

    c2 =  -0.824 + 4.717 / r_v
    c1 =  2.030 - 3.007 * c2

    x2 = x * x
    y = x2 - F99_X02
    d = x2 / (y * y + x2 * F99_GAMMA2)
    k = c1 + c2 * x + F99_C3 * d
    if x >= F99_C5:
        y = x - F99_C5
        y2 = y * y
        k += F99_C4 * (0.5392 * y2 + 0.05644 * y2 * y)

    return a_v * (1. + k / r_v)


_F99_XKNOTS = 1.e4 / np.array([np.inf, 26500., 12200., 6000., 5470.,
                               4670., 4110., 2700., 2600.])


def _f99_kknots(double[:] xknots, double r_v):
    cdef double c1, c2, d, x, x2, y, rv2
    cdef double[:] kknots_view
    cdef int i
    c2 =  -0.824 + 4.717 / r_v
    c1 =  2.030 - 3.007 * c2
    rv2 = r_v * r_v

    kknots = np.empty(9, dtype=np.float)
    kknots_view = kknots
    kknots_view[0] = -r_v
    kknots_view[1] = 0.26469 * r_v/3.1 - r_v
    kknots_view[2] = 0.82925 * r_v/3.1 - r_v
    kknots_view[3] = -0.422809 + 1.00270*r_v + 2.13572e-04*rv2 - r_v
    kknots_view[4] = -5.13540e-02 + 1.00216 * r_v - 7.35778e-05*rv2 - r_v
    kknots_view[5] = 0.700127 + 1.00184*r_v - 3.32598e-05*rv2 - r_v
    kknots_view[6] = (1.19456 + 1.01707*r_v - 5.46959e-03*rv2 +
                      7.97809e-04 * rv2 * r_v - 4.45636e-05 * rv2*rv2 - r_v)
    for i in range(7,9):
        x2 = xknots[i] * xknots[i]
        y = (x2 - F99_X02)
        d = x2 /(y * y + x2 * F99_GAMMA2)
        kknots_view[i] = c1 + c2*xknots[i] + F99_C3 * d

    return kknots


class Fitzpatrick99(object):
    """Fitzpatrick (1999) dust extinction function for arbitrary R_V.

    An instance of this class is a callable that can be used as
    ``f(wave, a_v)`` where ``wave`` is a 1-d array of wavelengths and ``a_v``
    is a scalar value.

    Parameters
    ----------
    r_v : float, optional
        R_V value. Default is 3.1.
    """

    def __init__(self, r_v=3.1):
        self.r_v = r_v

        kknots = _f99_kknots(_F99_XKNOTS, r_v)
        self._spline = splmake(_F99_XKNOTS, kknots, order=3)

    def __call__(self, np.ndarray wave not None, double a_v, unit='aa'):
        cdef double[:] wave_view, out_view
        cdef double r_v = self.r_v
        cdef double ebv = a_v / r_v
        cdef size_t i
        cdef size_t n

        # translate `wave` to inverse microns
        if unit == 'invum':
            pass
        elif unit == 'aa':
            wave = 1e4 / wave
        else:
            raise ValueError("unrecognized unit")

        # Optical/IR spline: evaluate at all wavelengths; we will overwrite
        # the UV points afterwards. These are outside the spline range, so
        # they don't actually take too much time to evaluate.
        out = spleval(self._spline, wave)  # this is actually "k"
        
        # Analytic function in the UV (< 2700 Angstroms).
        wave_view = wave
        out_view = out
        n = wave.shape[0]
        for i in range(n):
            # for optical/IR, `out` is actually k, but we wanted
            # a_v/r_v * (k+r_v), so we adjust here.
            if wave_view[i] < 1e4 / 2700.:
                out_view[i] = ebv * (out_view[i] + r_v)

            # for UV, we overwrite the array with the UV function value.
            else:
                out_view[i] = f99_uv_invum(wave_view[i], a_v, r_v)

        return out


# functional interface for Fitzpatrick (1999) with R_V = 3.1
_fitzpatrick99_fixed = F99(3.1)


def fitzpatrick99(wave, a_v, unit='aa'):
    """Fitzpatrick (1999) dust extinction law for R_V = 3.1.

    Parameters
    ----------
    wave : numpy.ndarray (1-d)
        Input wavelengths or wavenumbers (see units).
    a_v : float
        Total V-band extinction in magnitudes.
    unit : {'aa', 'invum'}, optional
        Wavelength units: Angstroms or inverse microns.

    Returns
    -------
    Extinction in magnitudes at each input wavelength.
    """
    return _fitzpatrick99_fixed(wave, a_v, unit=unit)


# -----------------------------------------------------------------------------
# Fitzpatrick & Massa 2007

DEF FM07_X0 = 4.592
DEF FM07_GAMMA = 0.922
DEF FM07_C1 = -0.175
DEF FM07_C2 = 0.807
DEF FM07_C3 = 2.991
DEF FM07_C4 = 0.319
DEF FM07_C5 = 6.097
DEF FM07_X02 = FM07_X0 * FM07_X0
DEF FM07_GAMMA2 = FM07_GAMMA * FM07_GAMMA
DEF FM07_R_V = 3.1  # Fixed for the time being (used in fm07kknots)

# Used for wave < 2700.
cdef inline double fm07_uv_invum(double x, double a_v):
    """fm07 for x < 1e4/2700 microns^-1 (ultraviolet)"""
    cdef double d, x2, y, k

    x2 = x * x
    y = x2 - FM07_X02
    d = x2 / (y*y + x2 * FM07_GAMMA2)
    k = FM07_C1 + FM07_C2 * x + FM07_C3 * d
    if x > FM07_C5:
        y = x - FM07_C5
        k += FM07_C4 * y * y

    return a_v * (1. + k / FM07_R_V)


def _fm07_kknots(double[:] xknots):
    cdef double d
    cdef int i, n

    n = xknots.shape[0]
    kknots = np.empty(n, dtype=np.float)
    for i in range(0, 5):
        kknots[i] = (-0.83 + 0.63*FM07_R_V) * xknots[i]**1.84 - FM07_R_V
    kknots[5] = 0.
    kknots[6] = 1.322
    kknots[7] = 2.055
    for i in range(8, 10):
        d = xknots[i]**2 / ((xknots[i]**2 - FM07_X02)**2 +
                            xknots[i]**2 * FM07_GAMMA2)
        kknots[i] = FM07_C1 + FM07_C2 * xknots[i] + FM07_C3 * d

    return kknots


_fm07_xknots = np.array([0., 0.25, 0.50, 0.75, 1., 1.e4/5530., 1.e4/4000.,
                        1.e4/3300., 1.e4/2700., 1.e4/2600.])
_fm07_spline = splmake(_fm07_xknots, _fm07_kknots(_fm07_xknots), order=3)


def fm07(double[:] wave, double a_v, unit='aa'):
    """Fitzpatrick & Massa (2007) extinction model for R_V = 3.1.

    The Fitzpatrick & Massa (2007) [1]_ model, which has a slightly
    different functional form from that of Fitzpatrick (1999) [3]_
    (`extinction_f99`). Fitzpatrick & Massa (2007) claim it is
    preferable, although it is unclear if signficantly so (Gordon et
    al. 2009 [2]_). Defined from 910 A to 6 microns.

    {0}


    References
    ----------
    .. [1] Fitzpatrick, E. L. & Massa, D. 2007, ApJ, 663, 320
    .. [2] Gordon, K. D., Cartledge, S., & Clayton, G. C. 2009, ApJ, 705, 1320
    .. [3] Fitzpatrick, E. L. 1999, PASP, 111, 63
    """

    cdef double[:] wave_view, out_view
    cdef double ebv = a_v / FM07_R_V
    cdef size_t i
    cdef size_t n

    # translate `wave` to inverse microns
    if unit == 'invum':
        pass
    elif unit == 'aa':
        wave = 1e4 / wave
    else:
        raise ValueError("unrecognized unit")

    # Optical/IR spline: evaluate at all wavelengths; we will overwrite
    # the UV points afterwards. These are outside the spline range, so
    # they don't actually take too much time to evaluate.
    out = spleval(_fm07_spline, wave)  # this is actually "k"
        
    # Analytic function in the UV (< 2700 Angstroms).
    wave_view = wave
    out_view = out
    n = wave.shape[0]
    for i in range(n):
        # for optical/IR, `out` is actually k, but we wanted
        # a_v/r_v * (k+r_v), so we adjust here.
        if wave_view[i] < 1e4 / 2700.:
            out_view[i] = ebv * (out_view[i] + FM07_R_V)

        # for UV, we overwrite the array with the UV function value.
        else:
            out_view[i] = fm07_uv_invum(wave_view[i], a_v)

    return out


# -----------------------------------------------------------------------------
# Calzetti 2000
# http://adsabs.harvard.edu/abs/2000ApJ...533..682C

cdef inline double calzetti00k_uv_invum(double x):
    """calzetti00 `k` for 0.12 microns < wave < 0.63 microns (UV/optical),
    x in microns^-1"""
    return 2.659 * (((0.011*x - 0.198)*x + 1.509)*x - 2.156)


cdef inline double calzetti00k_ir_invum(double x):
    """calzetti00 `k` for 0.63 microns < wave < 2.2 microns (optical/IR),
    x in microns^-1"""
    return 2.659 * (1.040*x - 1.857)


cdef inline double calzetti00_invum(double x, double r_v):
    cdef double k
    
    if x > 1.5873015873015872:  # 1. / 0.63
        k = calzetti00k_uv_invum(x)
    else:
        k = calzetti00k_ir_invum(x)

    return 1.0 + k / r_v


def calzetti00(double[:] wave, double a_v, double r_v, unit='aa',
               np.ndarray out=None):
    """calzetti00(wave, a_v, r_v, unit='aa', out=None)

    Calzetti (2000) extinction function.

    Calzetti et al. (2000, ApJ 533, 682) developed a recipe for
    dereddening the spectra of galaxies where massive stars dominate the
    radiation output, valid between 0.12 to 2.2 microns. They estimate
    R_V = 4.05 +/- 0.80 from optical-IR observations of 4 starburst galaxies.


    Parameters
    ----------
    wave : numpy.ndarray (1-d)
        Wavelengths or wavenumbers.
    a_v : float
        Scaling parameter, A_V: extinction in magnitudes at characteristic
        V band  wavelength.
    r_v : float
        Ratio of total to selective extinction, A_V / E(B-V).
    unit : {'aa', 'invum'}, optional
        Unit of wave: 'aa' (Angstroms) or 'invum' (inverse microns).
    out : np.ndarray, optional
        If specified, store output values in this array.

    Returns
    -------
    Extinction in magnitudes at each input wavelength.
    """

    cdef:
        size_t i
        size_t n = wave.shape[0]
        double a = 0.0
        double b = 0.0

    if out is None:
        out = np.empty(n, dtype=np.float)
    else:
        assert out.shape == wave.shape
        assert out.dtype == np.float

    cdef scalar_func convert_wave
    if unit == 'aa':
        convert_wave = &aa_to_invum
    elif unit == 'invum':
        convert_wave = &noop
    else:
        raise ValueError("unrecognized unit")

    cdef double[:] out_view = out
    for i in range(n):
        out_view[i] = a_v * calzetti00_invum(convert_wave(wave[i]), r_v)

    return out


# ------------------------------------------------------------------------------
# convenience function for applying extinction to flux values, optionally
# in-place. (It turns out that this isn't really faster than just doing it
# in pure python...)

def apply(double[:] extinction, np.ndarray flux not None, bint inplace=False):
    """apply(extinction, flux, inplace=False)

    Apply extinction in magnitudes to flux values, optionally in-place.

    The output value is flux * 10**(-0.4 * extinction).
    """

    cdef size_t i
    cdef size_t n = extinction.shape[0]
    cdef double[:] out_view
    cdef double[:] flux_view = flux

    if not extinction.shape[0] == flux.shape[0]:
        raise ValueError("shape of flux and extinction must match")

    if inplace:
        out = flux
    else:
        out = np.empty(n, dtype=np.float)

    out_view = out
    for i in range(n):
        out_view[i] = 10.**(-0.4 * extinction[i]) * flux_view[i]

    return out
