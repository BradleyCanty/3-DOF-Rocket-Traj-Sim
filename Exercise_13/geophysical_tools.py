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
OMEGA_EARTH = 7.29211510E-5 #rad/s, Earth's angular velocity along z-axis of ECI frame
r_earth = 6371E3 #m, mean radius of earth

#For U.S. Standard Atmosphere Calculations:
R_star = 8314.32 #m/(kmol K), universal gas constant
r_0 = 6356766 #m, effective radius of the earth
gamma = 1.400 #specific heat ratio
T_0 = 288.15 #K, sstandard temperature at sea level
g_0prime = 9.80665 #m^2/(s^2*m^prime), geopotential constant 
g_0 = 9.80665 #m^2/(s^2), gravitational-field strength at sea level
Gamma = 1 #m^prime/m, unit-converting constant
M_0 = 28.9644 #kg/kmol, mean molecular weight of air
P_0 = 101325.0 #Pa, standard pressure at sea level

'''
Atmospheric Properties
'''
#Table 4 from U.S. Standard Atmosphere 1976
#Column 0: H_b
#Column 1: L_M,b
#Column 2: T_M,b
#Column 3: P_b
table4=np.array([[00000 , -0.0065 , 288.150 , 1.01325000000000E+5],
                 [11000 , 0.0000 , 216.650 , 2.26320639734629E+4],
                 [20000 , 0.0010 , 216.650 , 5.47488866967777E+3],
                 [32000 , 0.0028 , 228.650 , 8.68018684755228E+2],
                 [47000 , 0.0000 , 270.650 , 1.10906305554966E+2],
                 [51000 , -0.0028 , 270.650 , 6.69388731186873E+1],
                 [71000 , -0.0020 , 214.650 , 3.95642042804073E+0],
                 [84852 , 0.0000 , 186.946 , 3.73383589976215E-1]])

def get_air_temp(h):
    """
    Given height in meters, returns temperature in Kelvin
    Temperature is valid for -5,000 m < h < 86,000 m

    Parameters
    ----------
    h : float
        Altitude in meters

    Returns
    -------
    float
        air temperature in Kelvin
    """
    if h < 86000:
        H = r_0*h/(r_0+h)
        for b in range(table4.shape[0]-1):
            if H < table4[b+1][0]:
                break
        return table4[b][2] + table4[b][1]*(H-table4[b][0])
    else: #h is above 86
        #Convert h to km instead of meters
        h /= 1000
        if h < 91:
            return 186.8673
        elif h < 110:
            return 263.1905 - 76.3232 * math.sqrt(1 - ((h - 91) / -19.9429)**2)
        elif h < 120:
            return 240 + 12 * (h - 110)
        elif h <= 1000:
            eta = (h - 120) * (6356.766 + 120) / (6356.766 + h)
            return 1000 - 640 * math.exp(-0.01875 * eta)
        else:
            eta = (h - 120) * (6356.766 + 120) / (6356.766 + 1000)
            return 1000 - 640 * math.exp(-0.01875 * eta)
    
def get_air_pressure(h):
    """
    Given altitude in meters, returns pressure in Pascals.
    Pressure is valid for -5,000 m < z < 86,000 m.

    Parameters
    ----------
    h : float
        altitude in meters

    Returns
    -------
    float
        pressure in Pascals
    """
    if h < 86000:
        H = h*Gamma*r_0/(r_0+h)
        for b in range(table4.shape[0]-1):
            if (H < table4[b+1][0]):
                break
        C = -g_0*M_0/R_star
        Hb = table4[b][0]
        Lb = table4[b][1]
        Tb = table4[b][2]
        Pb = table4[b][3]
        if abs(Lb)>1E-12:
            return Pb * math.pow(1+Lb/Tb*(H-Hb),C/Lb)
        else:
            return Pb * math.exp(C*(H-Hb)/Tb)
    else: #h is above 86
        #Convert h to km instead of meters
        h /= 1000
        if h < 91:
            A = 0.000000
            B = 2.159582E-06
            C = -4.836957E-04
            D = -0.1425192
            E = 13.47530
            return math.exp(A*h**4+B*h**3+C*h**2+D*h+E)
        elif h < 100:
            A = 0.000000
            B = 3.304895E-05
            C = -0.009062730
            D = 0.6516698
            E = -11.03037
            return math.exp(A*h**4+B*h**3+C*h**2+D*h+E)
        elif h <110:
            A = 0.000000
            B = 6.693926E-05
            C = -0.01945388
            D = 1.719080
            E = -47.75030
            return math.exp(A*h**4+B*h**3+C*h**2+D*h+E)
        elif h < 120:
            A = 0.000000
            B = -6.539316E-05
            C = 0.02485568
            D = -3.223620
            E = 135.9355
            return math.exp(A*h**4+B*h**3+C*h**2+D*h+E)
        elif h < 150:
            A = 2.283506E-07
            B = -1.343221E-04
            C = 0.02999016
            D = -3.055446
            E = 113.5764
            return math.exp(A*h**4+B*h**3+C*h**2+D*h+E)
        elif h < 200:
            A = 1.209434E-08
            B = -9.692458E-06
            C = 0.003002041
            D = -0.4523015
            E = 19.19151
            return math.exp(A*h**4+B*h**3+C*h**2+D*h+E)
        elif h < 300:
            A = 8.113942E-10
            B = -9.822568E-07
            C = 4.687616E-04
            D = -0.1231710
            E = 3.067409
            return math.exp(A*h**4+B*h**3+C*h**2+D*h+E)
        elif h < 500:
            A = 9.814674E-11
            B = -1.654439E-07
            C = 1.148115E-04
            D = -0.05431334	
            E = -2.011365
            return math.exp(A*h**4+B*h**3+C*h**2+D*h+E)
        elif h < 750:
            A = -7.835161E-11
            B = 1.964589E-07
            C = -1.657213E-04
            D = 0.04305869
            E = -14.77132
            return math.exp(A*h**4+B*h**3+C*h**2+D*h+E)
        elif h <= 1000:
            A = 2.813255E-11
            B = -1.120689E-07
            C = 1.695568E-04
            D = -0.1188941
            E = 14.56718
            return math.exp(A*h**4+B*h**3+C*h**2+D*h+E)
        else:
            return 0
        

def get_air_density(h):
    """
    SUMMARY.

    Parameters
    ----------
    h : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.
    """
    if h < 86000:
        T = get_air_temp(h)
        P = get_air_pressure(h)
        return P*M_0/(R_star*T)
    
    else: #h is above 86
        #Convert h to km instead of meters
        h /= 1000
        if h < 91:
            A = 0.000000
            B = -3.322622E-06
            C = 9.111460E-04
            D = -0.2609971
            E = 5.944694
            return math.exp(A*h**4+B*h**3+C*h**2+D*h+E)
        elif h < 100:
            A = 0.000000
            B = 2.873405E-05
            C = -0.008492037
            D = 0.6541179
            E = -23.62010
            return math.exp(A*h**4+B*h**3+C*h**2+D*h+E)
        elif h < 110:
            A = -1.240774E-05
            B = 0.005162063
            C = -0.8048342
            D = 55.55996
            E = -1443.338
            return math.exp(A*h**4+B*h**3+C*h**2+D*h+E)
        elif h < 120:
            A = 0.00000
            B = -8.854164E-05
            C = 0.03373254
            D = -4.390837
            E = 176.5294
            return math.exp(A*h**4+B*h**3+C*h**2+D*h+E)
        elif h < 150:
            A = 3.661771E-07
            B = -2.154344E-04
            C = 0.04809214
            D = -4.884744
            E = 172.3597
            return math.exp(A*h**4+B*h**3+C*h**2+D*h+E)
        elif h < 200:
            A = 1.906032E-08
            B = -1.527799E-05
            C = 0.004724294
            D = -0.6992340
            E = 20.50921
            return math.exp(A*h**4+B*h**3+C*h**2+D*h+E)
        elif h < 300:
            A = 1.199282E-09
            B = -1.451051E-06
            C = 6.910474E-04
            D = -0.1736220
            E = -5.321644
            return math.exp(A*h**4+B*h**3+C*h**2+D*h+E)
        elif h < 500:
            A = 1.140564E-10
            B = -2.130756E-07
            C = 1.570762E-04
            D = -0.07029296
            E = -12.89844
            return math.exp(A*h**4+B*h**3+C*h**2+D*h+E)
        elif h < 750:
            A = 8.105631E-12
            B = -2.358417E-09
            C = -2.635110E-06
            D = -0.01562608
            E = -20.02246
            return math.exp(A*h**4+B*h**3+C*h**2+D*h+E)
        elif h <= 1000:
            A = -3.701195E-12
            B = -8.608611E-09
            C = 5.118829E-05
            D = -0.06600998
            E = -6.137674
            return math.exp(A*h**4+B*h**3+C*h**2+D*h+E)
        else:
            return 0

def get_speed_of_sound(h):
    """
    SUMMARY.

    Parameters
    ----------
    h : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.
    """
    
    T = get_air_temp(h)
    return math.sqrt(gamma*R_star*T/M_0)

def get_mach_number(v_norm,h):
    """
    SUMMARY.

    Parameters
    ----------
    v_norm : TYPE
        DESCRIPTION.
    h : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.
    """
    
    a = get_speed_of_sound(h)
    return v_norm/a

if __name__ == '__main__':
    h = 0
    print(f'For h = {h} [m]...')
    print(f'T = {get_air_temp(h)}')
    print(f'P = {get_air_pressure(h)}')
    print(f'rho = {get_air_density(h)}')
    print(f'a = {get_speed_of_sound(h)}')
    
    h = 20000
    print(f'For h = {h} [m]...')
    print(f'T = {get_air_temp(h)}')
    print(f'P = {get_air_pressure(h)}')
    print(f'rho = {get_air_density(h)}')
    print(f'a = {get_speed_of_sound(h)}')
    
    h = 85500
    print(f'For h = {h} [m]...')
    print(f'T = {get_air_temp(h)}')
    print(f'P = {get_air_pressure(h)}')
    print(f'rho = {get_air_density(h)}')
    print(f'a = {get_speed_of_sound(h)}')
    
    h = 110000
    print(f'For h = {h} [m]...')
    print(f'T = {get_air_temp(h)}')
    print(f'P = {get_air_pressure(h)}')
    print(f'rho = {get_air_density(h)}')
    print(f'a = {get_speed_of_sound(h)}')
    
    h = 1000000
    print(f'For h = {h} [m]...')
    print(f'T = {get_air_temp(h)}')
    print(f'P = {get_air_pressure(h)}')
    print(f'rho = {get_air_density(h)}')
    print(f'a = {get_speed_of_sound(h)}')
    
    h = 1000001
    print(f'For h = {h} [m]...')
    print(f'T = {get_air_temp(h)}')
    print(f'P = {get_air_pressure(h)}')
    print(f'rho = {get_air_density(h)}')
    print(f'a = {get_speed_of_sound(h)}')
