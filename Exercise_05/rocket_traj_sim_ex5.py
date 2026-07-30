# -*- coding: utf-8 -*-
"""
rocket_traj_sim_ex5.py

Description:
    Compute air density as a function of altitude using the isothermal barometric
    formula. Compare this model to tabulated U.S. Standard Atmosphere data
    (https://www.engineeringtoolbox.com/standard-atmosphere-d_604.html#gsc.tab=0)
    to verify that it's suitable for our use.
    
    The isothermal barometric formula is given by:
    rho(h) = rho0 * exp(-h/H)
    where
    rho0 = 1.225 [kg/m^3] = air density at sea level
    h = altitude above sea level
    H = scale height of the atmosphere = 8500 [m]
    
"""
import numpy as np
from matplotlib import pyplot as plt
import os

img_dir = os.getcwd() + r'\Images'
#If the Image directory does not exist, create it
if not os.path.exists(img_dir):
    os.makedirs(img_dir)

#Form the isothermal barometric formula
rho0 = 1.225 #sea level air density
H = 8500 #scale height
rho_air_isoth = lambda h: rho0 * np.exp(-h/H) #kg/m^3

#U.S. Standard Atmosphere altitude data from 0 to 80[km]
h_data = np.array([0,
1000,
2000,
3000,
4000,
5000,
6000,
7000,
8000,
9000,
10000,
15000,
20000,
25000,
30000,
40000,
50000,
60000,
70000,
80000
])

#U.S. Standard Atmosphere air density data corresponding to altitudes from 0 to 80[km]
rho_air_data = np.array([1.225,
1.112,
1.007,
0.9093,
0.8194,
0.7364,
0.6601,
0.59,
0.5258,
0.4671,
0.4135,
0.1948,
0.08891,
0.04008,
0.01841,
0.003996,
0.001027,
0.0003097,
0.00008283,
0.00001846
])

#Compare the isothermal barometric formula with the U.S. Standard Atmosphere data
h = np.linspace(0,80000,101)
plt.figure()
plt.plot(rho_air_isoth(h),h,'b',label='isothermal barometric formula')
plt.plot(rho_air_data,h_data,'r.',label='U.S. Std. Atmos. data')
plt.title(r'Air Density vs Altitude')
plt.xlabel(r'air density [$kg/m^3$]')
plt.ylabel('altitude [m]')
plt.grid(True)
plt.legend()
plt.savefig(img_dir + '/air_density_vs_altitude.png', dpi = 600)
plt.show()

print(f"Plots output to:\n{img_dir}")












