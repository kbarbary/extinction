#!/usr/bin/env python
"""Plot extinction functions for comparison"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.axes_grid1 import make_axes_locatable

import extinction

rcParams['font.family'] = 'serif'

wave = np.logspace(np.log10(910.), np.log10(30000.), 2000)

a_lambda = {'ccm89': extinction.ccm89(wave, 1.0, 3.1),
            'odonnell94': extinction.odonnell94(wave, 1.0, 3.1),
            'fitzpatrick99': extinction.fitzpatrick99(wave, 1.0),
            'fm07': extinction.fm07(wave, 1.0)}

names = list(a_lambda.keys())  # consistent ordering between panels

fig = plt.figure(figsize=(8.5, 6.))

ax = plt.axes()
for name in names:
    plt.plot(wave, a_lambda[name], label=name)
plt.axvline(x=2700., ls=':', c='k')
plt.axvline(x=3030.3030, ls=':', c='k')
plt.axvline(x=9090.9091, ls=':', c='k')
plt.axvspan(wave[0], 1150., fc='0.8', ec='none', zorder=-1000)
plt.axvspan(1150., 1250., fc='0.9', ec='none', zorder=-1000)    
plt.text(0.5, 0.95, '$R_V = 3.1$', transform=ax.transAxes, va='top',
         size='x-large')
plt.ylabel('Extinction ($A(\lambda)$ / $A_V$)')
plt.legend()
plt.setp(ax.get_xticklabels(), visible=False)

divider = make_axes_locatable(ax)
axresid = divider.append_axes("bottom", size=2.0, pad=0.2, sharex=ax)
for name in names:
    plt.plot(wave, a_lambda[name] - a_lambda['fitzpatrick99'])
plt.axvline(x=2700., ls=':', c='k')
plt.axvline(x=3030.3030, ls=':', c='k')
plt.axvline(x=9090.9091, ls=':', c='k')
plt.axvspan(wave[0], 1150., fc='0.8', ec='none', zorder=-1000)
plt.axvspan(1150., 1250., fc='0.9', ec='none', zorder=-1000)
plt.xlim(wave[0], wave[-1])
plt.ylim(ymin=-0.1,ymax=0.4)
plt.ylabel('residual from fitzpatrick99')
plt.xlabel(r'Wavelength ($\mathrm{\AA}$)')

ax.set_xscale('log')
axresid.set_xscale('log')
plt.tight_layout()

fig.savefig("_static/extinction.png")
