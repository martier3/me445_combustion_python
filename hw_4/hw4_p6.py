#!/usr/bin/env python

"""
A freely-propagating, premixed methane flat flame with multicomponent
transport properties.

Requires: cantera >= 2.5.0
Keywords: combustion, 1D flow, premixed flame, multicomponent transport,
          saving output
"""

import cantera as ct
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Simulation parameters
p = ct.one_atm * 10  # pressure [Pa]
Tin = 300.0  # unburned gas temperature [K]
reactants = 'CH4:1, O2:2, N2:7.52'  # premixed gas composition
width = 0.03  # m
loglevel = 1  # amount of diagnostic output (0 to 8)

# Solution object used to compute mixture properties, set to the state of the
# upstream fuel-air mixture
gas = ct.Solution('gri30.xml')
gas.TPX = Tin, p, reactants

# Set up flame object
f = ct.FreeFlame(gas, width=width)
f.set_refine_criteria(ratio=3, slope=0.06, curve=0.12)
f.show_solution()

# Solve with mixture-averaged transport model
f.transport_model = 'Mix'
f.solve(loglevel=loglevel, auto=True)

# Solve with the energy equation enabled
try:
    # save to HDF container file if h5py is installed
    f.write_hdf('adiabatic_flame.h5', group='mix', mode='w',
                description='solution with mixture-averaged transport')
except ImportError:
    f.save('adiabatic_flame.yaml', 'mix',
           'solution with mixture-averaged transport')

f.show_solution()
print('mixture-averaged flamespeed = {0:7f} m/s'.format(f.velocity[0]))

# Solve with multi-component transport properties
f.transport_model = 'Multi'
f.solve(loglevel)  # don't use 'auto' on subsequent solves
f.show_solution()
print('multicomponent flamespeed = {0:7f} m/s'.format(f.velocity[0]))
try:
    f.write_hdf('adiabatic_flame.h5', group='multi',
                description='solution with multicomponent transport')
except ImportError:
    f.save('adiabatic_flame.yaml', 'multi',
           'solution with multicomponent transport')

# write the velocity, temperature, density, and mole fractions to a CSV file
f.write_csv('concentrations.csv', quiet=False)


csv_filename = 'concentrations.csv'
plot_filename = 'concentration_plot'
title = 'Methane-Air product species concentration vs temp'

data = pd.read_csv(csv_filename, delimiter=',', header=0)
data_array = data.to_numpy()

CH4_X = data_array[:, 17]
CO_X = data_array[:, 18]
CO2_X = data_array[:, 19]
H2O_X = data_array[:, 9]
O2_X = data_array[:, 7]

temperature = data_array[:, 2]

plt.figure(1)
plt.plot(temperature, CH4_X, linewidth=1, label='CH4_X')
plt.plot(temperature, CO_X, linewidth=1, label='CO_X')
plt.plot(temperature, CO2_X, linewidth=1, label='CO2_X')
plt.plot(temperature, H2O_X, linewidth=1, label='H2O_X')
plt.plot(temperature, O2_X, linewidth=1, label='O2_X')

plt.xlabel('Temperature (C)')
plt.ylabel('Concentration')
plt.legend()
plt.title(title)
plt.savefig(plot_filename)
