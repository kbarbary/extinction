extinction
==========

*Fast interstellar dust extinction laws in Python*

Cython-optimized implementations of empirical dust exitinction laws
found in the literature.

Installation
------------

Using conda::

    conda install -c conda-forge extinction

Using pip (requires a C compiler)::

    pip install extinction

Extinction depends on numpy and scipy.

**Development version / source code:** http://github.com/kbarbary/extinction


Usage
-----

**Functions:**

.. autosummary::
   :toctree: api

   extinction.ccm89
   extinction.odonnell94
   extinction.calzetti00
   extinction.fitzpatrick99
   extinction.fm07
   extinction.apply
   extinction.remove

**Classes:**

.. autosummary::
   :toctree: api

   extinction.Fitzpatrick99

**Examples:**

Get extinction in magnitudes at a set of wavelengths for various dust laws::

  >>> import numpy as np
  >>> import extinction

  >>> wave = np.array([2000., 4000., 8000.])  # wavelength in Angstroms

  # Cardelli, Clayton & Mathis (1989) with A_V = 1 and R_V = 3.1
  >>> extinction.ccm89(wave, 1.0, 3.1)
  array([ 2.84252644,  1.4645557 ,  0.59748901])  # extinction in magnitudes
   
  # O'Donnell (1994)
  >>> extinction.odonnell94(wave, 1.0, 3.1)
  array([ 2.84252644,  1.42617802,  0.60793495])

  # Fitzpatrick (1999)
  >>> extinction.fitzpatrick99(wave, 1.0, 3.1)
  array([ 2.76225609,  1.79674653,  1.42325373])


The Fitzpatrick & Massa (2007) function has a fixed :math:`R_V` of 3.1::

   >>> extinction.fm07(wave, 1.0)
   array([ 2.90478329,  1.42645161,  0.54703201])
   
All extinction laws accept a ``unit`` keyword to change the interpretation of
the wavelength array from Angstroms to inverse microns::

  >>> wave = np.array([5., 2.5, 1.25])  # wavelength in inverse microns

  >>> extinction.ccm89(wave, 1.0, 3.1, unit='invum')
  array([ 2.84252644,  1.4645557 ,  0.59748901])  # extinction in magnitudes


Redden or deredden
..................

To "redden" or "deredden" flux values by some amount, use the
``apply`` and ``remove`` convenience functions::


  >>> from extinction import ccm89, apply

  >>> flux = np.ones(3)

  # "redden" flux by A_V = 1.0
  >>> apply(ccm89(wave, 1.0, 3.1), flux)
  array([ 0.07294397,  0.25952412,  0.5767723 ])

  # "deredden" flux by A_V = 1.0
  >>> remove(ccm89(wave, 1.0, 3.1), flux)
  array([ 13.70915145,   3.85320647,   1.73378645])



Comparison of functions
.......................

.. plot::

   import numpy as np
   import extinction
   from extinction_plot import extinction_figure
   
   wave = np.logspace(np.log10(910.), np.log10(30000.), 2000)

   a_lambda = {'ccm89': extinction.ccm89(wave, 1.0, 3.1),
               'odonnell94': extinction.odonnell94(wave, 1.0, 3.1),
               'fitzpatrick99': extinction.fitzpatrick99(wave, 1.0),
               'fm07': extinction.fm07(wave, 1.0)}
   extinction_figure(wave, a_lambda, 'fitzpatrick99')

           
R_V dependence of odonnell94
............................

.. plot::

   from collections import OrderedDict
   import numpy as np
   import extinction
   from extinction_plot import extinction_figure
   
   wave = np.logspace(np.log10(910.), np.log10(30000.), 2000)

   a_lambda = OrderedDict([
       ('$R_V$ = 2.1', extinction.odonnell94(wave, 1.0, 2.1)),
       ('$R_V$ = 2.6', extinction.odonnell94(wave, 1.0, 2.6)),
       ('$R_V$ = 3.1', extinction.odonnell94(wave, 1.0, 3.1)),
       ('$R_V$ = 3.6', extinction.odonnell94(wave, 1.0, 3.6)),
       ('$R_V$ = 4.1', extinction.odonnell94(wave, 1.0, 4.1))
       ])
   extinction_figure(wave, a_lambda, '$R_V$ = 3.1',
                     residual_lims=(-0.2, 0.6),
                     title_text='odonnell94')



R_V dependence of Fitzpatrick99
...............................
.. plot::

   from collections import OrderedDict
   import numpy as np
   import extinction
   from extinction_plot import extinction_figure
   
   wave = np.logspace(np.log10(910.), np.log10(30000.), 2000)

   a_lambda = OrderedDict([
       ('$R_V$ = 2.1', extinction.Fitzpatrick99(2.1)(wave, 1.0)),
       ('$R_V$ = 2.6', extinction.Fitzpatrick99(2.6)(wave, 1.0)),
       ('$R_V$ = 3.1', extinction.Fitzpatrick99(3.1)(wave, 1.0)),
       ('$R_V$ = 3.6', extinction.Fitzpatrick99(3.6)(wave, 1.0)),
       ('$R_V$ = 4.1', extinction.Fitzpatrick99(4.1)(wave, 1.0))
       ])
   extinction_figure(wave, a_lambda, '$R_V$ = 3.1',
                     residual_lims=(-0.2, 0.6),
                     title_text='Fitzpatrick99')



A note on parameterization
..........................

Most extinction laws here have two parameters: :math:`A_V` and
:math:`R_V`.  :math:`A_V` is a simply a linear scaling parameter; that
is: ``ccm89(wave, 2.0, 3.1)`` is the same as ``2.0 * ccm89(wave, 1.0,
3.1)``. :math:`R_V` changes the *shape* of the extinction
curve, rather than just the amplitude.

Traditionally, the meaning ascribed to these parameters was that
:math:`A_V` is the extinction in the *V* band, and :math:`R_V`
describes the ratio of total to selective extinction: :math:`R_V = A_V
/ E(B-V)`, where :math:`E(B-V)` is the difference in extinction
between the *B* and *V* bands. While this is approximately correct,
the *measured* extinction of a source in actual *B* and *V* bandpasses
will generally depend on the source spectrum and the shape of the
specific bandpasses. Therefore, the :math:`A_V` and :math:`R_V` that
are parameters in our extinction law will not correspond perfectly to
measured *B* and *V* extinctions. So, in the context of these
extinction laws, it is best to think of :math:`A_V` and :math:`R_V` as
simply parameters that describe the amplitude and shape of the
wavelength dependence, rather than corresponding directly to measured
magnitudes.

Finally, for a given shape of the extinction curve (described by
:math:`R_V`), one can equally well use :math:`E(B-V)` as a linear
scaling parameter in place of :math:`A_V`, with the equivalence
:math:`E(B-V) = A_V / R_V`. Note that the above caution applies here:
:math:`E(B-V)` should be considered simply a parameter describing
amplitude of extinction; it will not correspond exactly to a
difference in measured *B* and *V* extinctions.


License and Credits
-------------------

The license is MIT. Part of this code originated in the specutils package.
