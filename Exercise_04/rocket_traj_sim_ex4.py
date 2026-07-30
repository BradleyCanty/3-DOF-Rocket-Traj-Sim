# -*- coding: utf-8 -*-
"""
rocket_traj_sim_ex4.py

Description:
    Using the results of exercise 3, add in the effects of drag. 
    Assume that a V2 has a diameter of 1.65m, a drag coefficient of 0.125 and 
    that the density of air is 1.22kg/m^3. Graph the resulting 1-D position and 
    velocity as a function of time. Calculate a maximum height and compare that 
    result to previous numbers.
    
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

tt,ry,vy,ay = get_rocket_state(ry0,vy0,t,'Euler')

plt.figure()
plt.plot(tt,ay,'r')
plt.title(r'ay vs Time for V2 Rocket (g = f(h), $\rho_{air} = 1.22\ [kg/m^3]$)')
plt.xlabel('time [s]')
plt.ylabel('acceleration [m/s^2]')
plt.grid(True)
plt.savefig(img_dir + '/V2_rocket_vs_time.png', dpi = 600)
plt.show()

plt.figure()
plt.plot(tt,vy,'r')
plt.title(r'Vy vs Time for V2 Rocket (g = f(h), $\rho_{air} = 1.22\ [kg/m^3]$)')
plt.xlabel('time [s]')
plt.ylabel('velocity [m/s]')
plt.grid(True)
plt.savefig(img_dir + '/V2_rocket_vs_time.png', dpi = 600)
plt.show()

plt.figure()
plt.plot(tt,ry-r_earth,'r')
plt.title(r'Altitude vs Time for V2 Rocket (g = f(h), $\rho_{air} = 1.22\ [kg/m^3]$)')
plt.xlabel('time [s]')
plt.ylabel('Altitude [m]]')
plt.grid(True)
plt.savefig(img_dir + '/V2_rocket_vs_time.png', dpi = 600)
plt.show()

print(f"Plots output to:\n{img_dir}")











