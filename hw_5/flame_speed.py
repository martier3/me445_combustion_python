'''
Instructions: 
- Move this file into a folder of your choice
- Review the variables and settings in the section with values to change
- Open an anaconda prompt (you may also use pycharm, spyder, or another ide of your choice)
- In the anaconda prompt type: ipython flame_speed.py
- The code will run and save the two output files containing the plots in the same folder as the code
- The default is to print all of the progress information if you just want to run it and not see 
   any information pop up until the end change the loglevel to 0
- This code will take a while to run. If it takes longer than 30-45 minutes please contact 
   the TA to obtain the resulting plots.
- If you have any questions feel free to contact the TA

This example is adapted from the flame_speed_with_sensitivity_analysis.ipynb
example provided on the cantera website. The example can be found here:
https://cantera.org/examples/jupyter/flames/flame_speed_with_sensitivity_analysis.ipynb.html
'''
# Import packages
import cantera as ct
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
# Set the plotting values for nice looking plots
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['legend.fontsize'] = 10
plt.rcParams['figure.figsize'] = (8,6)
# Get the best of both ggplot and seaborn
plt.style.use('ggplot')
plt.style.use('seaborn-colorblind')
plt.rcParams['figure.autolayout'] = True

###############################################################################
################## START OF THE SECTION WITH VALUES TO CHANGE #################
###############################################################################
#Inlet Temperature in Kelvin
To = 300

# Inlet Pressure in Pascals
Po = 101325 * 10

# Minimim value of phi to calculate
phi_min = 0.5

# Maximum value of phi to calculate
phi_max = 1.5

# Number of points to use for the equivalence ratios
number_of_points = 20

# File with the mechanism information
mechanism = 'gri30.cti'

# Chemical composition of the fuel
fuel = 'CH4'

# Chemical composition of the oxidizer
oxidizer = {'O2':1.0, 'N2':3.76}

# name of the file to save the flame speed plot in
flame_speed_file = 'flame_speed_equivalence_ratio_10atm.png'

# name of the file to save the sensitivity plot in
sensitivity_file = 'reaction_sensitivity_10atm.png'

# loglevel = 1 prints out solution progress
# loglevel = 0 does not print any progress information 
loglevel = 1
###############################################################################
################### END OF THE SECTION WITH VALUES TO CHANGE ##################
###############################################################################

# Create a dictionary to store flame speed values
flame_speeds = {}

setup_gas = ct.Solution(mechanism)
reaction_sensitivity = pd.DataFrame(index=setup_gas.reaction_equations(range(setup_gas.n_reactions)))
# Iterate over the indented code with equivalence ratio from phi_min
# to phi_max using the number_of_points evenly spaced in the range
phi_values = np.linspace(phi_min, phi_max, num=number_of_points)
# round the values
phi_values = [round(x, 2) for x in phi_values]
#Iterate and solve for each phi value
for phi in phi_values:
        # Create a cantera gas object with specified mechanism
        gas = ct.Solution(mechanism) 

        # Create the premixed  fuel/air mixture 
        gas.set_equivalence_ratio(phi, fuel, oxidizer)
        # Set the temperature and pressure
        gas.TP = To, Po

        # Domain width in metres
        width = 0.014

        # Create the flame object
        flame = ct.FreeFlame(gas, width=width)

        # Define tolerances for the solver
        flame.set_refine_criteria(ratio=3, slope=0.1, curve=0.1)

        # Sove the flame with the default settings. 
        flame.solve(loglevel=loglevel, auto=True)
        # Obtain the flame speed. 
        # The speed used is the inlet velocity of the gas entering the domain
        Su0 = flame.velocity[0]
        # Print out the flame speed, equivalence ratio, and max temperature
        print("Flame Speed is: {:.2f} cm/s".format(Su0*100),
              'phi={:.2f}'.format(phi),
              'T_max={:.2f}'.format(max(flame.T)))
        # Save the flame speed to the corresponding spot in the dictionary
        flame_speeds[phi] = Su0
        reaction_sensitivity[phi] = flame.get_flame_speed_reaction_sensitivities()

# Create a matplotlib figure
plt.figure()
# Plot the flame speed as a function of equivalence ratio with a solid line
plt.plot(flame_speeds.keys(), flame_speeds.values(), '-')
# Add axis labels
plt.xlabel(r'$\Phi$')
plt.ylabel('Flame Speed (m/s)')
# Save the figure
plt.savefig(flame_speed_file, dpi=300)


plt.figure(figsize=(6, 12))
# Get in sensitivity values for the n most sensitive reactions
sensitive_reactions = reaction_sensitivity.nlargest(15,
                                                    reaction_sensitivity.columns,
                                                    keep='all')
# Plot the values for every other phi value to reduce clutter                                     
sensitive_reactions[sensitive_reactions.columns[::2]].plot.barh()
# Add the xlabel
plt.xlabel(r'Sensitivity: $\frac{\partial\:\ln{S_{u}}}{\partial\:\ln{k}}$')
# Save the figure, Higher dpi is used here since there are quite a few bars 
plt.savefig(sensitivity_file, dpi=500)
