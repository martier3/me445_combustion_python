"""
Adiabatic flame temperature and equilibrium composition for a fuel/air mixture
as a function of equivalence ratio, including formation of solid carbon.

Requires: cantera >= 2.5.0, matplotlib >= 2.0
Keywords: equilibrium, combustion, multiphase
"""

import matplotlib.pyplot as plt
import cantera as ct
import numpy as np
import csv

##############################################################################
# Edit these parameters to change the initial temperature, the pressure, and
# the phases in the mixture.

T = 300.0
P = 101325.0 * 10

# phases
gas = ct.Solution('gri30.yaml')

# the phases that will be included in the calculation, and their initial moles
mix_phases = [(gas, 1.0)]

# gaseous fuel species
fuel_species = 'CH4'

# equivalence ratio range
npoints = 500
phi = np.linspace(0.3, 3.5, npoints)

##############################################################################

mix = ct.Mixture(mix_phases)

# create some arrays to hold the data
tad = np.zeros(npoints)
xeq = np.zeros((mix.n_species, npoints))

for i in range(npoints):
    # set the gas state
    gas.set_equivalence_ratio(phi[i], fuel_species, 'O2:2.0, N2:7.52')

    # create a mixture of 1 mole of gas
    mix = ct.Mixture(mix_phases)
    mix.T = T
    mix.P = P

    # equilibrate the mixture adiabatically at constant P
    mix.equilibrate('HP', solver='gibbs', max_steps=1000)

    tad[i] = mix.T
    print('At phi = {0:12.4g}, Tad = {1:12.4g}'.format(phi[i], tad[i]))
    xeq[:, i] = mix.species_moles

# write output CSV file for importing into Excel
csv_file = 'adiabatic.csv'
with open(csv_file, 'w', newline='') as outfile:
    writer = csv.writer(outfile)
    writer.writerow(['phi', 'T (K)'] + mix.species_names)
    for i in range(npoints):
        writer.writerow([phi[i], tad[i]] + list(xeq[:, i]))
print('Output written to {0}'.format(csv_file))

plt.figure(1)
plt.plot(phi, tad)
plt.xlabel('Equivalence Ratio')
plt.ylabel('Temperature (K)')
plt.title('Equivalence Ratio vs. Adiabatic Temp (K)')
plt.savefig('adiabatic_temp.png')
