#!python
#cython: boundscheck=False, wraparound=False, initializedcheck=False, cdivision=True

"""Interstellar dust extinction functions."""

import numpy as np
cimport numpy as np

__version__ = "0.4.2"

__all__ = ['ccm89', 'odonnell94', 'Fitzpatrick99', 'fitzpatrick99', 'fm07',
           'calzetti00', 'apply', 'remove']


# We use some C code for Cubic splines in the Fitzpatrick99 and fm07 functions.
# There are two reasons for using this over something from scipy:
#
# First, Fitzpatrick99 specifies "natural" boundary conditions on the splines,
# but in scipy.interpolate, this can only be achieved using the CubicSpline
# class, which is new in scipy v0.18 (a bit too recent).
#
# Second, this C code is significantly faster than CubicSpline.
include "extern/bsplines.pxi"

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

    The parameters given in the original paper [1]_ are used.
    The claimed validity is 1250 Angstroms to 3.3 microns.

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
    In Cardelli, Clayton & Mathis (1989) the mean
    R_V-dependent extinction law, is parameterized as

    .. math::

       <A(\lambda)/A_V> = a(x) + b(x) / R_V

    where the coefficients a(x) and b(x) are functions of
    wavelength. At a wavelength of approximately 5494.5 angstroms (a
    characteristic wavelength for the V band), a(x) = 1 and b(x) = 0,
    so that A(5494.5 angstroms) = A_V. This function returns

    .. math::

       A(\lambda) = A_V (a(x) + b(x) / R_V)

    References
    ----------
    .. [1] Cardelli, J. A., Clayton, G. C., & Mathis, J. S. 1989, ApJ, 345, 245

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
    """odonnell94(wave, a_v, r_v, unit='aa', out=None)

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


cdef double *_F99_XKNOTS = [0.0, 1.e4 / 26500., 1.e4 / 12200., 1.e4 / 6000.,
                            1.e4 / 5470., 1.e4 / 4670., 1.e4 / 4110.,
                            1.e4 / 2700., 1.e4 / 2600.]


cdef void _f99_kknots(double *xknots, double r_v, double *out):
    """Fill knot k values in Fitzpatrick 99 for given R_V.
    xknots and out should be arrays of length 9."""

    cdef double c1, c2, d, x, x2, y, rv2
    cdef int i

    c2 =  -0.824 + 4.717 / r_v
    c1 =  2.030 - 3.007 * c2
    rv2 = r_v * r_v

    out[0] = -r_v
    out[1] = 0.26469 * r_v/3.1 - r_v
    out[2] = 0.82925 * r_v/3.1 - r_v
    out[3] = -0.422809 + 1.00270*r_v + 2.13572e-04*rv2 - r_v
    out[4] = -5.13540e-02 + 1.00216 * r_v - 7.35778e-05*rv2 - r_v
    out[5] = 0.700127 + 1.00184*r_v - 3.32598e-05*rv2 - r_v
    out[6] = (1.19456 + 1.01707*r_v - 5.46959e-03*rv2 +
              7.97809e-04 * rv2 * r_v - 4.45636e-05 * rv2*rv2 - r_v)
    for i in range(7,9):
        x2 = xknots[i] * xknots[i]
        y = (x2 - F99_X02)
        d = x2 / (y * y + x2 * F99_GAMMA2)
        out[i] = c1 + c2*xknots[i] + F99_C3 * d


cdef class Fitzpatrick99(object):
    """Fitzpatrick (1999) dust extinction function for arbitrary R_V.

    An instance of this class is a callable that can be used as
    ``f(wave, a_v)`` where ``wave`` is a 1-d array of wavelengths and ``a_v``
    is a scalar value.

    Parameters
    ----------
    r_v : float, optional
        R_V value. Default is 3.1.

    Examples
    --------
    Create a callable that gives the extinction law for a given ``r_v``
    and use it:

    >>> f = Fitzpatrick99(3.1)

    >>> f(np.array([3000., 4000.]), 1.)
    array([ 1.79939955,  1.42338583])

    """

    cdef bs_spline1d *spline   # pointer to c struct
    cdef readonly double r_v
    
    def __cinit__(self, double r_v=3.1):
        self.r_v = r_v

        cdef double *kknots = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        _f99_kknots(_F99_XKNOTS, r_v, kknots)

        cdef bs_bcs bcs
        cdef bs_exts exts
        cdef bs_errorcode code

        bcs.left.type = BS_DERIV2
        bcs.left.value = 0.0
        bcs.right.type = BS_DERIV2
        bcs.right.value = 0.0
        exts.left.type = BS_VALUE
        exts.left.value = 0.0
        exts.right.type = BS_VALUE
        exts.right.value = 0.0
        code = bs_spline1d_create(bs_array(_F99_XKNOTS, 9, 1),
                                  bs_array(kknots, 9, 1),
                                  bcs, exts, &self.spline)
        assert_ok(code)


    def __call__(self, np.ndarray wave not None, double a_v, unit='aa'):
        cdef double[:] wave_view, out_view
        cdef double r_v = self.r_v
        cdef double ebv = a_v / r_v
        cdef size_t i
        cdef size_t n
        cdef bs_errorcode code
        
        # translate `wave` to inverse microns
        if unit == 'invum':
            pass
        elif unit == 'aa':
            wave = 1e4 / wave
        else:
            raise ValueError("unrecognized unit")

        n = wave.shape[0]
        out = np.empty(n, dtype=np.float64)

        wave_view = wave
        out_view = out

        # Optical/IR spline: evaluate at all wavelengths; we will overwrite
        # the UV points afterwards. These are outside the spline range, so
        # they don't actually take too much time to evaluate.
        # (outview actually holds "k" after this call; see below)
        code = bs_spline1d_eval(self.spline, to_bs_array(wave_view),
                                to_bs_array(out_view))
        assert_ok(code)
        
        # Analytic function in the UV (< 2700 Angstroms).
        for i in range(n):
            # for optical/IR, `out` is actually k, but we wanted
            # a_v/r_v * (k+r_v), so we adjust here.
            if wave_view[i] < 3.7037037037037037:  # 1e4 / 2700
                out_view[i] = ebv * (out_view[i] + r_v)

            # for UV, we overwrite the array with the UV function value.
            else:
                out_view[i] = f99_uv_invum(wave_view[i], a_v, r_v)

        return out

    def __dealloc__(self):
        bs_spline1d_free(self.spline)


# functional interface for Fitzpatrick (1999) with R_V = 3.1
_fitzpatrick99_fixed = Fitzpatrick99(3.1)


def fitzpatrick99(wave, a_v, r_v=3.1, unit='aa'):
    """fitzpatrick99(wave, a_v, r_v=3.1, unit='aa')

    Fitzpatrick (1999) dust extinction function.

    Fitzpatrick (1999) [1]_ model which relies on the parametrization
    of Fitzpatrick & Massa (1990) [2]_ in the UV (below 2700 A) and
    spline fitting in the optical and IR. This function is defined
    from 910 A to 6 microns, but note the claimed validity goes down
    only to 1150 A. The optical spline points are not taken from F99
    Table 4, but rather updated versions from E. Fitzpatrick (this
    matches the Goddard IDL astrolib routine FM_UNRED).


    Parameters
    ----------
    wave : numpy.ndarray (1-d)
        Input wavelengths or wavenumbers (see units).
    a_v : float
        Total V-band extinction in magnitudes.
    r_v : float
        Ratio of total to selective extinction, A_V / E(B-V).
    unit : {'aa', 'invum'}, optional
        Wavelength units: Angstroms or inverse microns.

    Returns
    -------
    Extinction in magnitudes at each input wavelength.

    References
    ----------
    .. [1] Fitzpatrick, E. L. 1999, PASP, 111, 63
    .. [2] Fitzpatrick, E. L. & Massa, D. 1990, ApJS, 72, 163
    """
    if r_v == 3.1:
        return _fitzpatrick99_fixed(wave, a_v, unit=unit)
    else:
        return Fitzpatrick99(r_v)(wave, a_v, unit=unit)

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


cdef void _fm07_kknots(double *xknots, double *out):
    """Return k knots given xknots for fm07. (Could be a static array)"""
    cdef double d
    cdef int i

    for i in range(0, 5):
        out[i] = (-0.83 + 0.63*FM07_R_V) * xknots[i]**1.84 - FM07_R_V
    out[5] = 0.0
    out[6] = 1.322
    out[7] = 2.055
    for i in range(8, 10):
        d = xknots[i]**2 / ((xknots[i]**2 - FM07_X02)**2 +
                            xknots[i]**2 * FM07_GAMMA2)
        out[i] = FM07_C1 + FM07_C2 * xknots[i] + FM07_C3 * d


cdef double *_FM07_XKNOTS = [0., 0.25, 0.50, 0.75, 1., 1.e4 / 5530.,
                             1.e4 / 4000., 1.e4 / 3300., 1.e4 / 2700.,
                             1.e4 / 2600.]
cdef double *_FM07_KKNOTS = [0., 0., 0., 0., 0.,
                             0., 0., 0., 0., 0.]
_fm07_kknots(_FM07_XKNOTS, _FM07_KKNOTS)


def fm07(np.ndarray wave not None, double a_v, unit='aa'):
    """fm07(wave, a_v, unit='aa')

    Fitzpatrick & Massa (2007) extinction model for R_V = 3.1.

    The Fitzpatrick & Massa (2007) [1]_ model, which has a slightly
    different functional form from that of Fitzpatrick (1999) [3]_
    (`extinction_f99`). Fitzpatrick & Massa (2007) claim it is
    preferable, although it is unclear if signficantly so (Gordon et
    al. 2009 [2]_). Defined from 910 A to 6 microns.

    Parameters
    ----------
    wave : numpy.ndarray (1-d)
        Wavelengths or wavenumbers.
    a_v : float
        Scaling parameter, A_V: extinction in magnitudes at characteristic
        V band  wavelength.
    unit : {'aa', 'invum'}, optional
        Unit of wave: 'aa' (Angstroms) or 'invum' (inverse microns).

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

    # create the spline (we do this here rather than globally so we can
    # deallocate it.)
    cdef bs_spline1d* spline
    cdef bs_bcs bcs
    cdef bs_exts exts
    cdef bs_errorcode code

    bcs.left.type = BS_DERIV2
    bcs.left.value = 0.0
    bcs.right.type = BS_DERIV2
    bcs.right.value = 0.0
    exts.left.type = BS_VALUE
    exts.left.value = 0.0
    exts.right.type = BS_VALUE
    exts.right.value = 0.0
    code = bs_spline1d_create(bs_array(_FM07_XKNOTS, 10, 1),
                              bs_array(_FM07_KKNOTS, 10, 1),
                              bcs, exts, &spline)
    assert_ok(code)

    n = wave.shape[0]
    out = np.empty(n, dtype=np.float64)

    wave_view = wave
    out_view = out

    # Optical/IR spline: evaluate at all wavelengths; we will overwrite
    # the UV points afterwards. These are outside the spline range, so
    # they don't actually take too much time to evaluate.
    # (outview actually holds "k" after this call; see below)
    code = bs_spline1d_eval(spline, to_bs_array(wave_view),
                            to_bs_array(out_view))
    assert_ok(code)
    bs_spline1d_free(spline)
    
    # Analytic function in the UV (< 2700 Angstroms).
    for i in range(n):
        # for optical/IR, `out` is actually k, but we wanted
        # a_v/r_v * (k+r_v), so we adjust here.
        if wave_view[i] < 3.7037037037037037:  # 1e4 / 2700
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
    :math:`R_V = 4.05 \pm 0.80` from optical-IR observations of
    4 starburst galaxies.

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
# convenience functions for applying/removing extinction to/from flux
# values, optionally in-place. It turns out that pure Python is about as
# fast as C here.

def apply(extinction, flux, inplace=False):
    """apply(extinction, flux, inplace=False)

    Apply extinction to flux values.

    This is a convenience function to perform "reddening" of
    flux values. It simply performs ``flux * 10**(-0.4 * extinction)``:
    flux is decreased (for positive extinction values).

    Parameters
    ----------
    extinction : numpy.ndarray (1-d)
        Extinction in magnitudes. (Positive extinction values decrease flux
        values.)
    flux : numpy.ndarray
        Flux values.
    inplace : bool, optional
        Whether to perform the operation in-place on the flux array. If True,
        the return value is a reference to the input flux array.

    Returns
    -------
    new_flux : numpy.ndarray (1-d)
        Flux values with extinction applied.
    """

    trans = 10.**(-0.4 * extinction)

    if inplace:
        flux *= trans
        return flux
    else:
        return flux * trans

# ------------------------------------------------------------------------------
# convenience function for removing extinction from flux values, optionally
# in-place. (It turns out that this isn't really faster than just doing it
# in pure python...)

def remove(extinction, flux, inplace=False):
    """remove(extinction, flux, inplace=False)

    Remove extinction from flux values.

    This is a convenience function to "deredden" fluxes. It simply performs 
    ``flux * 10**(0.4 * extinction)``: flux is increased (for positive
    extinction values).

    Parameters
    ----------
    extinction : numpy.ndarray (1-d)
        Extinction in magnitudes.
    flux : numpy.ndarray
        Flux values.
    inplace : bool, optional
        Whether to perform the operation in-place on the flux array. If True,
        the return value is a reference to the input flux array.

    Returns
    -------
    new_flux : numpy.ndarray (1-d)
        Flux values with extinction removed.
    """

    trans = 10.**(0.4 * extinction)

    if inplace:
        flux *= trans
        return flux
    else:
        return flux * trans
