# -*- coding: utf-8 -*-
"""
rocket_traj_sim_ex1.py

Description:
    Formulate a 1D model for rocket motion in the absence of outside forces. Using 
    Euler method for numerical integration, compute the velocity of a V2 rocket and 
    compare it with the velocity computed using Tsiolkovsky's rocket equation.
    See 'Exercise_01_Notes.pdf' for solution steps.
    
"""
import numpy as np
import math
from matplotlib import pyplot as plt
import os

img_dir = os.getcwd() + r'\Images'
#If the Image directory does not exist, create it
if not os.path.exists(img_dir):
    os.makedirs(img_dir)
    
m_empty = 4000 #kg, empty mass
m_propellant = 8800 #kg
m_initial = m_empty + m_propellant
m_dot = 129.4 #kg/s, mass flow rate of fuel
v_exhaust = 2050 #m/s, exhaust velocity

#Initial conditions
t0 = 0
v0 = 0 
a0 = 0

'''
Acceleration given by 
a = T/m = v_exhaust * m_dot / (m_empty + m_fuel - m_dot * t)

Euler method of numerical integration given by:
v_i+1 = a_i * dt + v_i
'''
t = np.linspace(0, 10, 51)
v = np.zeros(len(t))
v[0] = v0
a = np.zeros(len(t))
a[0] = v_exhaust * m_dot / (m_empty + m_propellant)

#Assume uniform time step
dt = t[1] - t[0]

for ti in range (len(t)-1):
    a[ti+1] = v_exhaust * m_dot / (m_empty + m_propellant - m_dot*dt*(ti+1))
    v[ti+1] = a[ti]*dt + v[ti]
    #print(f't = {t[ti]}')
    
#Compute the theoretical deltaV at each time step
t2 = np.linspace(0,10,11)
m_final = lambda t: m_initial - m_dot * t
deltaV = np.zeros(len(t2))
deltaV[0] = v0
for i in range (len(t2)-1):
    deltaV[i+1] = v_exhaust * math.log(m_initial/m_final(i+1))

plt.figure()
plt.plot(t,v,'r',label='Euler method solution')
plt.plot(t2,deltaV,'b.',markersize=15,label='analytical solution')
plt.title('Velocity vs Time for V2 Rocket (no external forces)')
plt.xlabel('time [s]')
plt.ylabel('velocity [m/s]')
plt.legend()
plt.grid(True)
plt.savefig(img_dir + '/V2_vel_vs_time.png', dpi = 600)
plt.show()

print(f"Plots output to:\n{img_dir}")











