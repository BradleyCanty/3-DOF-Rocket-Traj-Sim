# -*- coding: utf-8 -*-
"""
rocket_traj_sim_ex6.py

Description:
    Fit speed of sound versus altitude data with piecewise continuous lines.
    Additionally, fit drag coefficient vs Mach number data (derived from V2 wind 
    tunnel tests) with a natural cubic spline.

    Mach number is defined as:
    M = v / a
    where
    a = speed of sound = sqrt(gamma*R*T)
    gamma = ratio of specific heats = cp/cv ~ 1.4
    R = specific gas constant of air = 287.05 [J/(kg*K)]
    T(h) = molecular temperature, which is a nonlinear function of altitude
         
    The goal here is to compute the speed of sound, a = sqrt(gamma*R*T)
    The specific heat ratio, gamma, can be assumed to be constant. 
    The specific gas constant of air, R, can be assumed to be constant.
    Therefore, whats left to be solved for is the molecular temperature, T.
    This can be modeled using piecewise linear temperature gradients:
    T(h) = T_i + L_i(h-H_i)
    where 
    i is the Layer index
    T_i = the temperature at the base of the atmosphere layer [K]
    L_i = the rate at which temperature changes with altitude (i.e. "lapse rate") [K/km]
    h = the altitude above Earth's surface (i.e. "geodetic altitude")
    H_i = the altitude at which the atmosphere layer begins (i.e. "base altitude") [km]
         
    | Layer | Base Altitude (H_i)| Base Temp (T_i) | Lapse Rate (L_i)   | Layer Description      |
    |       | [km]               | [K]             | [K/km]             |                        |
    +-------+--------------------+-----------------+--------------------+------------------------+
    | 0     |   0.0              | 288.15          | -6.5               | Troposphere            |
    | 1     |  11.0              | 216.65          |  0.0               | Tropopause             |
    | 2     |  20.0              | 216.65          |  1.0               | Stratosphere (Lower)   |
    | 3     |  32.0              | 228.65          |  2.8               | Stratosphere (Upper)   |
    | 4     |  47.0              | 270.65          |  0.0               | Stratopause            |
    | 5     |  51.0              | 270.65          | -2.8               | Mesosphere (Lower)     |
    | 6     |  71.0              | 214.65          | -2.0               | Mesosphere (Upper)     |
    | 7     |  84.852            | 186.87          |  0.0               | Mesopause (to 86 km)   |
    
    Write a function which computes the temperature at a geodetic height h. Use an if 
    statement to determine its atmosphere layer, then plug in the associated variables
    into the temperature equation and return the result
    
    Additionally, given Cd vs Mach number data for the V2 rocket [Ref. 1], fit the 
    data to a curve using a natural cubic spline [Ref. 2]. Then, plot both the curve 
    and the data.
    
    | Mach |  Cd |
    +------+----+
    | 0.0  | .25 |
    | 0.5  | .18 |
    | 0.7  | .17 |
    | 1.0  | .29 |
    | 1.15 | .38 |
    | 1.5  | .27 |
    | 1.75 | .21 |
    | 2.0  | .18 |
    | 2.5  | .15 |
    | 3.5  | .12 |
    | 5.0  | .1  |

Created by: Bradley Canty, 2026/06/21
    
References:
    1) Natural cubic spline algorithm: 
        'Numerical Analysis' by Burden and Faires, page 147
        
    2) Cd vs Re data for V2 rocket: 
        'A Ballistic Missile Primer' by Steve Fetter, page 7
        'A History of German Guided Missile Development' by Agardograph, page 63
"""
import os
import numpy as np
from matplotlib import pyplot as plt
from natural_cubic_spline_fxns import get_cubic_spline_points

img_dir = os.getcwd() + r'\Images'
#If the Image directory does not exist, create it
if not os.path.exists(img_dir):
    os.makedirs(img_dir)

#Function to get air temperature at altitude h
def get_temp(h):
    if h < 11000:      #Trophsphere
        return 288.15 - 6.5E-3 * (h - 0)
    elif h < 20000:    #Tropopause
        return 216.65
    elif h < 32000:    #Stratosphere (lower)
        return 216.65 + 1E-3 * (h - 20000)
    elif h < 47000:    #Stratosphere (upper)
        return 228.65 + 2.8E-3 * (h - 32000)
    elif h < 51000:     #Stratopause
        return 270.65 
    elif h < 71000:    #Mesophere (lower)
        return 270.65 - 2.8E-3 * (h - 51000)
    elif h < 84852:   #Mesophere (upper)
        return 214.65 - 2.0E-3 * (h - 71000)
    else:       #Mesopause (to 86 km)
        return 186.87

gamma = 1.4
R = 287.05 #J/(kg*K)
h = np.linspace(0,100000,1001) #m, altitude
T = np.zeros(len(h)) #K, temperature
a = np.zeros(len(h)) #m/s, speed of sound

#Compute air temperature and speed of sound corresponding to each altitude h
for i in range(len(h)):
    T[i] = get_temp(h[i])
    a[i] = np.sqrt(gamma*R*T[i])
    

M_data = [0.0,0.5,0.7,1,1.15,1.5,1.75,2.0,2.5,3.5,5]
Cd_data = [.25,.18,.17,.29,.38,.27,.21,.18,.15,.12,.1]
n = 101

#Generate the Cd vs Mach number profile
M_values,Cd_values = get_cubic_spline_points(M_data,Cd_data,n)

'''
GENERATE PLOTS
'''
plt.figure()
plt.plot(T,h,'r')
plt.title('Altitude vs Temperature')
plt.xlabel('temperature [K]')
plt.ylabel('altitude [m]')
plt.grid(True)
plt.savefig(img_dir + '/alt_vs_temp.png', dpi = 600)
plt.show()

plt.figure()
plt.plot(a,h,'r')
plt.title('Altitude vs Speed of Sound')
plt.xlabel('speed of sound [m/s]')
plt.ylabel('altitude [m]')
plt.grid(True)
plt.savefig(img_dir + '/alt_vs_speed_of_sound.png', dpi = 600)
plt.show()

plt.figure()
plt.plot(M_values,Cd_values,'r')
plt.title('Drag Coefficient vs Mach Number')
plt.xlabel('mach number')
plt.ylabel('drag coefficient')
plt.grid(True)
plt.savefig(img_dir + '/cd_vs_mach.png', dpi = 600)
plt.show()

print(f"Plots output to:\n{img_dir}")







