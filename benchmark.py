#!/usr/bin/env python
from __future__ import print_function, division

import sys
import timeit

args = sys.argv[1:]
if len(args) == 0:
    args = ['ccm89', 'fitzpatrick99', 'apply']

if 'ccm89' in args:
    setup = """
import numpy
from extinction import {name}
wave = numpy.logspace(numpy.log10({minwave}),
                      numpy.log10({maxwave}),
                      {size:d})
"""

    stmt = "{name}(wave, 1.0, 3.1)"

    benchmarks = [
        {'name': 'ccm89',
         'minwave': 1000., 'maxwave': 30000., 'size': 10, 'nloops': 1000000},
        {'name': 'ccm89',
         'minwave': 1000., 'maxwave': 30000., 'size': 100, 'nloops': 100000},
        {'name': 'ccm89',
         'minwave': 1000., 'maxwave': 30000., 'size': 1000, 'nloops': 10000},
        {'name': 'ccm89',
         'minwave': 1000., 'maxwave': 30000., 'size': 10000, 'nloops': 1000},
        {'name': 'odonnell94',
         'minwave': 1000., 'maxwave': 30000., 'size': 100, 'nloops': 100000},
        {'name': 'ccm89',
         'minwave': 1000., 'maxwave': 1e4/8, 'size': 100, 'nloops': 100000},
        {'name': 'ccm89',
         'minwave': 1e4/8, 'maxwave': 1e4/3.3, 'size': 100, 'nloops': 100000},
        {'name': 'ccm89',
         'minwave': 1e4/3.3, 'maxwave': 1e4/1.1, 'size': 100, 'nloops': 100000},
        {'name': 'ccm89',
         'minwave': 1e4/1.1, 'maxwave': 1e4/0.3, 'size': 100, 'nloops': 100000}]
    for b in benchmarks:
        times = timeit.repeat(stmt.format(**b), setup=setup.format(**b),
                              repeat=3, number=b['nloops'])
        b['time'] = min(times) / b['nloops'] * 1e6
        print('{name:5s}  wave=[{minwave: 7.1f}, {maxwave: 7.1f}]  size={size:5d}  {time:8.2f}us'.format(**b))


    print()

if 'fitzpatrick99' in args:
    
    # -----------------------------------------------------------------------
    # Fitzpatrick init
    setup = """
import numpy
from extinction import Fitzpatrick99
"""

    stmt = "Fitzpatrick99(3.1)"

    nloops = 10000
    times = timeit.repeat(stmt, setup=setup, repeat=3, number=nloops)
    t = min(times) / nloops * 1e6
    print('F99 init              {:8.2f}us'.format(t))


    #------------------------------------------------------------------------
    # Fitzpatrick fast method

    setup = """
import numpy
from extinction import Fitzpatrick99

wave = numpy.logspace(numpy.log10({minwave}),
                      numpy.log10({maxwave}),
                      {size:d})
f = Fitzpatrick99(3.1)
"""

    stmt = "f(wave, 1.0)"

    benchmarks = [
        {'minwave': 1000., 'maxwave': 30000., 'size': 10, 'nloops': 10000},
        {'minwave': 1000., 'maxwave': 30000., 'size': 100, 'nloops': 10000},
        {'minwave': 1000., 'maxwave': 30000., 'size': 1000, 'nloops': 10000},
        {'minwave': 1000., 'maxwave': 30000., 'size': 10000, 'nloops': 1000},
        {'minwave': 1000., 'maxwave': 2700.,  'size': 1000, 'nloops': 10000},
        {'minwave': 2700., 'maxwave': 30000., 'size': 1000, 'nloops': 10000}]
    for b in benchmarks:
        times = timeit.repeat(stmt.format(**b), setup=setup.format(**b),
                              repeat=3, number=b['nloops'])
        b['time'] = min(times) / b['nloops'] * 1e6
        print('Fitzpatrick99 call  wave=[{minwave: 7.1f}, {maxwave: 7.1f}]  size={size:5d}  {time:8.2f}us'.format(**b))


if 'apply' in args:

    setup = """
import numpy
from extinction import apply
ext = numpy.ones({size:d})
flux = numpy.ones({size:d})
"""

    stmt = "apply(ext, flux, inplace={inplace!r})"

    benchmarks = [
        {'size': 10, 'nloops': 100000, 'inplace': True},
        {'size': 10, 'nloops': 100000, 'inplace': False},
        {'size': 100, 'nloops': 100000, 'inplace': True},
        {'size': 100, 'nloops': 100000, 'inplace': False},
        {'size': 1000, 'nloops': 10000, 'inplace': False}]
    for b in benchmarks:
        times = timeit.repeat(stmt.format(**b), setup=setup.format(**b),
                              repeat=3, number=b['nloops'])
        b['time'] = min(times) / b['nloops'] * 1e6
        print('apply  size={size:5d} inplace={inplace!r} {time:8.2f}us'.format(**b))


#def benchmark(f, args, repeat=3):
#    """Time `f(*args)`, repeat `repeat` times and return the best time."""
#
#    number = 1
#    times = []
#
#    # determine number of times to run the loop.
#    while number < 10**9:
#        t0 = time.time()
#        for i in range(number):
#            f(*args)
#        t = time.time() - t0
#
#        if t > 0.1:
#            break
#
#        number *= 10
#
#    # run the loop `number` times
    
