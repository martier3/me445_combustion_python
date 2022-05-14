#!/usr/bin/env python

"""
Calculate steady-state solutions for a combustor, modeled as a single well-stirred
reactor, for different residence times.

We are interested in the steady-state burning solution. This example explores
the effect of changing the residence time on completeness of reaction (through
the burned gas temperature) and on the total heat release rate.

Demonstrates the use of a MassFlowController where the mass flow rate function
depends on variables other than time by capturing these variables from the
enclosing scope. Also shows the use of a PressureController to create a constant
pressure reactor with a fixed volume.

Requires: cantera >= 2.5.0, matplotlib >= 2.0
Keywords: combustion, reactor network, well-stirred reactor, plotting
"""

import numpy as np
import matplotlib.pyplot as plt
import cantera as ct

# Use reaction mechanism GRI-Mech 3.0. For 0-D simulations,
# no transport model is necessary.
gas = ct.Solution('gri30.yaml')

# Create a Reservoir for the inlet, set to a methane/air mixture at a specified
# equivalence ratio
equiv_ratio = 0.5  # lean combustion
gas.TP = 300.0, ct.one_atm
gas.set_equivalence_ratio(equiv_ratio, 'CH4:1.0', 'O2:1.0, N2:3.76')
inlet = ct.Reservoir(gas)

# Create the combustor, and fill it initially with a mixture consisting of the
# equilibrium products of the inlet mixture. This state corresponds to the state
# the reactor would reach with infinite residence time, and thus provides a good
# initial condition from which to reach a steady-state solution on the reacting
# branch.
gas.equilibrate('HP')
combustor = ct.IdealGasReactor(gas)
combustor.volume = 1.0

# Create a reservoir for the exhaust
exhaust = ct.Reservoir(gas)

# Use a variable mass flow rate to keep the residence time in the reactor
# constant (residence_time = mass / mass_flow_rate). The mass flow rate function
# can access variables defined in the calling scope, including state variables
# of the Reactor object (combustor) itself.


def mdot(t):
    return combustor.mass / residence_time


inlet_mfc = ct.MassFlowController(inlet, combustor, mdot=mdot)

# A PressureController has a baseline mass flow rate matching the 'master'
# MassFlowController, with an additional pressure-dependent term. By explicitly
# including the upstream mass flow rate, the pressure is kept constant without
# needing to use a large value for 'K', which can introduce undesired stiffness.
outlet_mfc = ct.PressureController(combustor, exhaust, master=inlet_mfc, K=0.01)

# the simulation only contains one reactor
sim = ct.ReactorNet([combustor])

# Run a loop over decreasing residence times, until the reactor is extinguished,
# saving the state after each iteration.
states = ct.SolutionArray(gas, extra=['tres'])

residence_time = 0.1  # starting residence time
while combustor.T > 300:
    sim.set_initial_time(0.0)  # reset the integrator
    sim.advance_to_steady_state()
    print('tres = {:.2e}; T = {:.1f}'.format(residence_time, combustor.T))
    states.append(combustor.thermo.state, tres=residence_time)
    residence_time *= 0.9  # decrease the residence time for the next iteration

CH4 = states.X[1:, 13]
O2 = states.X[1:, 3]
H2O = states.X[1:, 5]
CO2 = states.X[1:, 15]

# Plot results

temperature = states.T[:-1]

plt.figure(1)
plt.plot(temperature, CH4, linewidth=1, label='CH4')
plt.plot(temperature, O2, linewidth=1, label='O2')
plt.plot(temperature, H2O, linewidth=1, label='H2O')
plt.plot(temperature, CO2, linewidth=1, label='CO2')

plt.xlim([1300, 1500])

plt.xlabel('Temperature (C)')
plt.ylabel('Concentration')
plt.legend()
plt.title('Inlet Temperature vs Concentration')
plt.savefig('p7.png')
