# -*- coding: utf-8 -*-
"""
geophysical_tools.py
"""

import math
import numpy as np

'''
Physical Constants
'''
G = 6.674E-11 #m^3/(kg*s^2), universal gravitational constant
m_earth = 5.972E24 #kg, mass of earth
r_earth = 6371E3 #m, mean radius of earth
OMEGA_EARTH = 7.29211510E-5 #rad/s, Earth's angular velocity along z-axis of ECI frame
#OMEGA_EARTH = 0

'''
Atmospheric Properties
'''
def get_air_temp(h):
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

def get_speed_of_sound(h):
    gamma = 1.4
    R = 287.05 #J/(kg*K)
    T = get_air_temp(h)
    return math.sqrt(gamma*R*T)

def get_air_density(h):
    rho0 = 1.225 #sea level air density
    H = 8500 #scale height
    return rho0 * np.exp(-h/H) #kg/m^3

def get_mach_number(v_norm,h):
    a = get_speed_of_sound(h)
    return v_norm/a


    