# -*- coding: utf-8 -*-
"""
atmosphere_fxns.py

"""
import math
import numpy as np

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
