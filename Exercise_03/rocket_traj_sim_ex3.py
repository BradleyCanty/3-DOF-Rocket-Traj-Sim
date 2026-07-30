# -*- coding: utf-8 -*-
"""
rocket_traj_sim_ex3.py

Description:
    Modify the results from Exercise 2 to account for the weakening of gravity 
    as a function of distance from the center of the Earth. Graph the resulting 
    1-D position and velocity as a function of time. Additionally, solve the
    system of ODEs using the midpoint method (RK2) and the 4th Order Runge Kutta
    method.
    
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
plt.title('V2 Rocket $a_y$ vs Time (g = f(h))')
plt.xlabel('time [s]')
plt.ylabel('acceleration [m/s^2]')
plt.grid(True)
plt.savefig(img_dir + '/V2_ay_vs_time.png', dpi = 600)
plt.show()

plt.figure()
plt.plot(tt,vy,'r')
plt.title('V2 Rocket $v_y$ vs Time (g = f(h))')
plt.xlabel('time [s]')
plt.ylabel('velocity [m/s]')
plt.grid(True)
plt.savefig(img_dir + '/V2_vy_vs_time.png', dpi = 600)
plt.show()

plt.figure()
plt.plot(tt,ry-r_earth,'r')
plt.title('V2 Rocket Altitude vs Time(g = f(h))')
plt.xlabel('time [s]')
plt.ylabel('altitude [m]]')
plt.grid(True)
plt.savefig(img_dir + '/V2_rocket_vs_time.png', dpi = 600)
plt.show()

print(f"Plots output to:\n{img_dir}")











