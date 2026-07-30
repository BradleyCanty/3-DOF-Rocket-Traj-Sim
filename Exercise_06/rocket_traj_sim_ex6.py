# -*- coding: utf-8 -*-
"""
rocket_traj_sim_ex6.py

Description:
    Repeat Exercise 5 using the isothermal barometric formula for air density as
    a function of altitude. Plot the resulting 1-D acceleration, velocity, and position 
    as a function of time.
    
    The isothermal barometric formula is given by:
    rho(h) = rho0 * exp(-h/H)
    where
    rho0 = 1.225 [kg/m^3] = air density at sea level
    h = altitude above sea level
    H = scale height of the atmosphere = 8500 [m]
    
"""
import numpy as np
from matplotlib import pyplot as plt
from math_fxns import get_rocket_state
import os

img_dir = os.getcwd() + r'\Images'
#If the Image directory does not exist, create it
if not os.path.exists(img_dir):
    os.makedirs(img_dir)

r_earth = 6371000 #m, mean earth radius

#Initial conditions
ry0 = r_earth
vy0 = 0 

t = np.linspace(0, 500, 301)

#Assume uniform time step
dt = t[1] - t[0]

tt,ry,vy,ay = get_rocket_state(ry0,vy0,t,'RK4')

plt.figure()
plt.plot(tt,ay,'r')
plt.title(r'$a_y$ vs Time for V2 Rocket (g = f(h), $\rho_{air} = f(h)\ [kg/m^3]$)')
plt.xlabel('time [s]')
plt.ylabel('acceleration [m/s^2]')
plt.grid(True)
plt.savefig(img_dir + '/V2_rocket_vs_time.png', dpi = 600)
plt.show()

plt.figure()
plt.plot(tt,vy,'r')
plt.title(r'$v_y$ vs Time (g = f(h), $\rho_{air} = f(h)\ [kg/m^3]$)')
plt.xlabel('time [s]')
plt.ylabel('velocity [m/s]')
plt.grid(True)
plt.savefig(img_dir + '/V2_rocket_vs_time.png', dpi = 600)
plt.show()

plt.figure()
plt.plot(tt,ry-r_earth,'r')
plt.title(r'Altitude vs Time (g = f(h), $\rho_{air} = f(h)\ [kg/m^3]$)')
plt.xlabel('time [s]')
plt.ylabel('Altitude [m]]')
plt.grid(True)
plt.savefig(img_dir + '/V2_rocket_vs_time.png', dpi = 600)
plt.show()

print(f"Plots output to:\n{img_dir}")


