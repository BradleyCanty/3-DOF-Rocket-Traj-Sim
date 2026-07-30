 # -*- coding: utf-8 -*-
"""
rocket_traj_sim_ex8.py

Description:
    Implement the Mach number-varying drag coefficient in the 1D rocket trajectory
    simulation. Note that you must use the temperature vs altitude piecewise linear 
    curve fit from the previous exercise to determine the temperature as a function
    of altitude, and then use the temperature in the calculation of speed of sound.
    Velocity is computed at each time step in the simulation, while Mach number is
    the velocity divided by speed of sound, M = v/a. Once you have Mach number, use
    the Cd vs Mach cubic curve fit to find the the corresponding drag coefficient.

    Plot the resulting 1-D position, velocity, and acceleration as a function of time.
    Additionally, plot the Mach number and drag coefficient as a function of time.

Created by: Bradley Canty, 2026/06/21
"""
import os
import numpy as np
from matplotlib import pyplot as plt
from natural_cubic_spline_fxns import get_cubic_spline_point
from atmosphere_fxns import get_speed_of_sound

img_dir = os.getcwd() + r'\Images'
#If the Image directory does not exist, create it
if not os.path.exists(img_dir):
    os.makedirs(img_dir)

m_empty = 4000 #kg, empty mass
m_propellant = 8800 #kg
m_initial = m_empty + m_propellant
m_dot = 129.4 #kg/s, mass flow rate of fuel
v_exhaust = 2050 #m/s, exhaust velocity
t_burnout = m_propellant/m_dot #seconds
g0 = 9.81
r_earth = 6371000 #m, mean earth radius
g = lambda h: g0*r_earth**2/(r_earth + h)**2
rho0 = 1.225 #sea level air density
H = 8500 #scale height
rho_air = lambda h: rho0 * np.exp(-h/H) #kg/m^3
D = 1.65 #m, diameter of V2 rocket
A = np.pi*(D/2)**2

def get_mach_number(v,h):
    a = get_speed_of_sound(h)
    return np.abs(v)/a

#Data for Cd vs Mach profile
M_data = [0.0,0.5,0.7,1,1.15,1.5,1.75,2.0,2.5,3.5,5]
Cd_data = [.25,.18,.17,.29,.38,.27,.21,.18,.15,.12,.1]
def get_Cd(M):
    if M > M_data[-1]:
        #If the Mach number is greater than 5, assume the Cd is equivalent to that at Mach 5
        return Cd_data[-1]
    else: 
        return get_cubic_spline_point(M_data,Cd_data,M)

#Initial conditions
t0 = 0
a0 = 0
v0 = 0 
y0 = r_earth

'''
Acceleration given by 
a = T/m = v_exhaust * m_dot / (m_empty + m_fuel - m_dot * t)

Euler-Cromer method of numerical integration given by:
v_i+1 = a_i+1 * dt + v_i
'''
t = np.linspace(0, 500, 601)

#Record the components of acceleration for plotting

h = np.zeros(len(t))
h[0] = y0 - r_earth

a_grav = np.zeros(len(t))
a_grav[0] = -g(h[0])
a_thrust = np.zeros(len(t))
a_thrust[0] = v_exhaust * m_dot / (m_empty + m_propellant)
a_drag = np.zeros(len(t))
a_drag[0] = 0

ay = np.zeros(len(t))
ay[0] = a_grav[0] + a_thrust[0] + a_drag[0]
vy = np.zeros(len(t))
vy[0] = v0
ry = np.ones(len(t)) * r_earth
ry[0] = y0


M = np.zeros(len(t))
Cd = np.zeros(len(t))

#Assume uniform time step
dt = t[1] - t[0]

for i in range (len(t)-1):
    a_grav[i+1] = - g(ry[i] - r_earth)
    M[i] = get_mach_number(vy[i], h[i])
    Cd[i] = get_Cd(M[i])
    if (t[i] < t_burnout):
        a_thrust[i+1] = v_exhaust * m_dot / (m_empty + m_propellant - m_dot*dt*(i+1))
        a_drag[i+1] = - (1/2 * rho_air(h[i]) * vy[i] * abs(vy[i]) * A * Cd[i]) / (m_empty + m_propellant - m_dot*dt*(i+1))
        
    else:
        a_thrust[i+1] = 0
        a_drag[i+1] = - (1/2 * rho_air(h[i]) * vy[i] * abs(vy[i]) * A * Cd[i]) / m_empty

    ay[i+1] = a_grav[i+1] + a_thrust[i+1] + a_drag[i+1]
    vy[i+1] = ay[i+1]*dt + vy[i]
    ry[i+1] = vy[i+1]*dt + ry[i]
    h[i+1] = ry[i+1] - r_earth
    
    if h[i+1] <= 0: break

'''
GENERATE PLOTS
'''
plt.figure()
plt.plot(t,ay,'r')
plt.title(r'$a_y$ vs Time for V2 Rocket (g = f(h), $\rho_{air} = f(h) [kg/m^3]$)')
plt.xlabel('time [s]')
plt.ylabel('acceleration [m/s^2]')
plt.grid(True)
plt.savefig(img_dir + '/V2_ay_vs_time.png', dpi = 600)
plt.show()

plt.figure()
plt.plot(t,a_thrust,'r')
plt.title(r'$a_y,thrust$ vs Time for V2 Rocket (g = f(h), $\rho_{air} = f(h) [kg/m^3]$)')
plt.xlabel('time [s]')
plt.ylabel('acceleration [m/s^2]')
plt.grid(True)
plt.savefig(img_dir + '/V2_ay_thrust_vs_time.png', dpi = 600)
plt.show()

plt.figure()
plt.plot(t,vy,'r')
plt.title(r'$v_y$ vs Time (g = f(h), $\rho_{air} = f(h) [kg/m^3]$)')
plt.xlabel('time [s]')
plt.ylabel('velocity [m/s]')
plt.grid(True)
plt.savefig(img_dir + '/V2_vy_vs_time.png', dpi = 600)
plt.show()

plt.figure()
plt.plot(t,h,'r')
plt.title(r'Altitude vs Time (g = f(h), $\rho_{air} = f(h) [kg/m^3]$)')
plt.xlabel('time [s]')
plt.ylabel('Altitude [m]')
plt.grid(True)
plt.savefig(img_dir + '/V2_alt_vs_time.png', dpi = 600)
plt.show()

plt.figure()
plt.plot(t,a_drag,'r')
plt.title(r'Drag Acceleration vs Time (g = f(h), $\rho_{air} = f(h) [kg/m^3]$)')
plt.xlabel('time [s]')
plt.ylabel(r'drag acceleration [$m/s^2$]')
plt.grid(True)
plt.savefig(img_dir + '/V2_drag_accel_vs_time.png', dpi = 600)
plt.show()

plt.figure()
plt.plot(t, M,'r',label='Mach number')
plt.plot(t,Cd,'b',label='drag coefficient')
plt.title('Mach Number and Drag Coefficient vs Time')
plt.legend()
plt.grid(True)
plt.savefig(img_dir + '/V2_mach_and_cd_vs_time.png', dpi = 600)
plt.show()

print(f"Plots output to:\n{img_dir}")

