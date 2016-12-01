"""Plot extinction functions for comparison"""

from collections import OrderedDict

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from mpl_toolkits.axes_grid1 import make_axes_locatable

import extinction

rcParams['font.family'] = 'serif'


def extinction_figure(wave, a_lambda, residual_from, residual_lims=(-0.1, 0.4), title_text='$R_V = 3.1$'):

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
    plt.text(0.65, 0.95, title_text, transform=ax.transAxes, va='top',
             ha='right', size='x-large')
    plt.ylabel('Extinction ($A(\lambda)$ / $A_V$)')
    plt.legend()
    plt.setp(ax.get_xticklabels(), visible=False)

    divider = make_axes_locatable(ax)
    axresid = divider.append_axes("bottom", size=2.0, pad=0.2, sharex=ax)
    for name in names:
        plt.plot(wave, a_lambda[name] - a_lambda[residual_from])
    plt.axvline(x=2700., ls=':', c='k')
    plt.axvline(x=3030.3030, ls=':', c='k')
    plt.axvline(x=9090.9091, ls=':', c='k')
    plt.axvspan(wave[0], 1150., fc='0.8', ec='none', zorder=-1000)
    plt.axvspan(1150., 1250., fc='0.9', ec='none', zorder=-1000)
    plt.xlim(wave[0], wave[-1])
    plt.ylim(ymin=residual_lims[0], ymax=residual_lims[1])
    plt.ylabel('residual from ' + residual_from)
    plt.xlabel(r'Wavelength ($\mathrm{\AA}$)')

    ax.set_xscale('log')
    axresid.set_xscale('log')
    plt.tight_layout()

    return fig
