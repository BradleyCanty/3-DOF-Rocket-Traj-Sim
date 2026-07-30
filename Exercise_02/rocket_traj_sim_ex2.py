# -*- coding: utf-8 -*-
"""
rocket_traj_sim_ex2.py

Description:
    Modify the model from Exercise 1 by adding in a constant gravity. Note that 
    while g will not change with time, the weight will. Write code so that when 
    the rocket runs out of fuel, its mass stops changing, and the thrust goes to 
    zero. Additionally, put the Euler method into a separate 'math_fxns.py' file
    and import it into this script. Graph the resulting 1-D position, velocity 
    and acceleration as a function of time.
    
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
g0 = 9.81 #m/s^2
rho_air = 1.225 #kg/m^3

#Initial conditions
ry0 = r_earth
vy0 = 0 

m_propellant = 8800 #kg
m_dot = 129.4 #kg/s, mass flow rate of fuel
v_exhaust = 2050 #m/s, exhaust velocity
t_burnout = m_propellant/m_dot #seconds

#Specify times to evaluate the trajectory
t1 = np.linspace(0, t_burnout, 301)
t2 = np.linspace(t_burnout, 300, 301)

tt1,ry1,vy1,ay1,m1 = get_rocket_state(ry0,vy0,t1,'Euler')
tt2,ry2,vy2,ay2,m2 = get_rocket_state(ry1[-1],vy1[-1],t2,'Euler')

#Post processing
ay1[-1] = ay2[0]
m1[-1] = m2[0]
ay2[-1] = ay2[-2]
m2[-1] = m2[-2]

#calculate dynamic pressure
q1 = rho_air*vy1**2/2
q2 = rho_air*vy2**2/2

plt.figure()
plt.plot(tt1,q1,'r',label='main engine burn')
plt.plot(tt2,q2,'b',label='coast')
plt.title('V2 Dynamic Pressure vs Time (g = 9.81[$m/s^2$])')
plt.xlabel('time [s]')
plt.ylabel(r'dynamic pressure [$N/m^2$]')
plt.legend()
plt.grid(True)
plt.savefig(img_dir + '/V2_dynamic_pres_vs_time.png', dpi = 600)
plt.show()

plt.figure()
plt.plot(tt1,(ry1-r_earth)/1000,'r',label='main engine burn')
plt.plot(tt2,(ry2-r_earth)/1000,'b',label='coast')
plt.title('V2 Altitude vs Time (g = 9.81[$m/s^2$])')
plt.xlabel('time [s]')
plt.ylabel('altitude [km]]')
plt.legend()
plt.grid(True)
plt.savefig(img_dir + '/V2_altitude_vs_time.png', dpi = 600)
plt.show()

plt.figure()
plt.plot(tt1,ay1/g0,'r',label='main engine burn')
plt.plot(tt2,ay2/g0,'b',label='coast')
plt.title('V2 Acceleration (gs) vs Time (g = 9.81[$m/s^2$])')
plt.xlabel('time [s]')
plt.ylabel(r'acceleration [gs]')
plt.legend()
plt.grid(True)
plt.savefig(img_dir + '/V2_accel_vs_time.png', dpi = 600)
plt.show()

plt.figure()
plt.plot(tt1,vy1,'r',label='main engine burn')
plt.plot(tt2,vy2,'b',label='coast')
plt.title('V2 $V_y$ vs Time (g = 9.81[$m/s^2$])')
plt.xlabel('time [s]')
plt.ylabel('velocity [m/s]')
plt.legend()
plt.grid(True)
plt.savefig(img_dir + '/V2_vy_vs_time.png', dpi = 600)
plt.show()

plt.figure()
plt.plot(tt1,m1,'r',label='main engine burn')
plt.plot(tt2,m2,'b',label='coast')
plt.title('V2 Mass vs Time (g = 9.81[$m/s^2$])')
plt.xlabel('time [s]')
plt.ylabel('mass [kg]')
plt.legend()
plt.grid(True)
plt.savefig(img_dir + '/V2_mass_vs_time.png', dpi = 600)
plt.show()

print(f"Plots output to:\n{img_dir}")











