# -*- coding: utf-8 -*-
"""
geophysical_tools.py
"""

import math
import numpy as np
import coordinate_conversion_tools as cct

'''
-------------------------------------------------------------------------------
PHYSICAL CONSTANTS
-------------------------------------------------------------------------------
'''
#EARTH MODEL (WGS84)
G = 6.674E-11 #m^3/(kg*s^2), universal gravitational constant
m_earth = 5.972E24 #kg, mass of earth
OMEGA_EARTH = 7.2921151467e-5 #rad/s, Earth's angular velocity along z-axis of ECI frame
R_eq = 6378137.0 # km, Earth equatorial radius
f = 1 / 298.257223563 # flattening
R_p = R_eq*(1-f) # km, Earth polar radius
e2 = 2*f - f**2 #first eccentricity squared
eprime2 = (R_eq**2 - R_p**2)/R_p**2 #second eccentricity squared
#R_mean = (2*R_eq + R_p)/3 #m, mean radius of Earth

#ATMOSPHERE MODEL (U.S. STANDARD ATMOSPHERE 1976)
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
-------------------------------------------------------------------------------
ATMOSPHERE FUNCTIONS
-------------------------------------------------------------------------------
'''
#Table 4 from U.S. Standard Atmosphere 1976
#Column           0: H_b  1: L_M,b  2: T_M,b  3: P_b
table4=np.array([[00000 , -0.0065,  288.150,  1.01325000000000E+5],
                 [11000 , 0.0000 ,  216.650,  2.26320639734629E+4],
                 [20000 , 0.0010 ,  216.650,  5.47488866967777E+3],
                 [32000 , 0.0028 ,  228.650,  8.68018684755228E+2],
                 [47000 , 0.0000,   270.650,  1.10906305554966E+2],
                 [51000 , -0.0028,  270.650,  6.69388731186873E+1],
                 [71000 , -0.0020,  214.650,  3.95642042804073E+0],
                 [84852 , 0.0000,   186.946,  3.73383589976215E-1]])

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

'''
-------------------------------------------------------------------------------
GEODETIC FUNCTIONS
-------------------------------------------------------------------------------
'''
def vincenty_inverse(lat1: float, lon1: float, lat2: float, lon2: float, 
                     max_iter: int = 200, tol: float = 1e-12) -> dict:
    """
    Calculates geodesic distance and azimuths between two points on WGS-84 ellipsoid.
    
    Parameters:
      lat1, lon1: Latitude/Longitude of Point 1 in degrees
      lat2, lon2: Latitude/Longitude of Point 2 in degrees
      max_iter  : Maximum convergence iterations
      tol       : Convergence tolerance for lambda (longitude difference on sphere)
      
    Returns:
      dict containing:
        - 'distance_m': Geodesic distance in meters
        - 'azimuth_initial_deg': Forward azimuth at start (0-360)
        - 'azimuth_final_deg': Forward azimuth at end (0-360)
        - 'converged': Boolean flag
    """
    # 1. WGS-84 Ellipsoid Constants
    a = 6378137.0          # Semi-major axis (meters)
    f = 1 / 298.257223563   # Flattening
    b = (1 - f) * a        # Semi-minor axis (meters)

    # Coincident points check
    if abs(lat1 - lat2) < 1e-9 and abs(lon1 - lon2) < 1e-9:
        return {'distance_m': 0.0, 'azimuth_initial_deg': 0.0, 
                'azimuth_final_deg': 0.0, 'converged': True}

    # Convert degrees to radians
    phi1, L1 = math.radians(lat1), math.radians(lon1)
    phi2, L2 = math.radians(lat2), math.radians(lon2)

    # Reduced latitudes (u1, u2)
    U1 = math.atan((1 - f) * math.tan(phi1))
    U2 = math.atan((1 - f) * math.tan(phi2))
    sinU1, cosU1 = math.sin(U1), math.cos(U1)
    sinU2, cosU2 = math.sin(U2), math.cos(U2)

    # Difference in longitude
    L = L2 - L1
    lambda_lon = L  # Initial estimate for lambda

    # Iterative loop
    for iteration in range(max_iter):
        sin_lambda = math.sin(lambda_lon)
        cos_lambda = math.cos(lambda_lon)

        # Angular distance sigma on auxiliary sphere
        sin_sigma = math.sqrt(
            (cosU2 * sin_lambda)**2 + 
            (cosU1 * sinU2 - sinU1 * cosU2 * cos_lambda)**2
        )
        
        # Coincident or antipodal edge cases
        if sin_sigma == 0:
            return {'distance_m': 0.0, 'azimuth_initial_deg': 0.0, 
                    'azimuth_final_deg': 0.0, 'converged': True}

        cos_sigma = sinU1 * sinU2 + cosU1 * cosU2 * cos_lambda
        sigma = math.atan2(sin_sigma, cos_sigma)

        # Forward azimuth alpha of the geodesic at the equator
        sin_alpha = (cosU1 * cosU2 * sin_lambda) / sin_sigma
        cos2_alpha = 1 - sin_alpha**2

        # Angular distance sigma_m from equator to geodesic midpoint
        if cos2_alpha == 0:
            cos2_sigma_m = 0.0  # Equatorial line
        else:
            cos2_sigma_m = cos_sigma - (2 * sinU1 * sinU2) / cos2_alpha

        # Correction C
        C = (f / 16) * cos2_alpha * (4 + f * (4 - 3 * cos2_alpha))

        # Update lambda
        lambda_prev = lambda_lon
        lambda_lon = L + (1 - C) * f * sin_alpha * (
            sigma + C * sin_sigma * (
                cos2_sigma_m + C * cos_sigma * (-1 + 2 * cos2_sigma_m**2)
            )
        )

        # Check convergence
        if abs(lambda_lon - lambda_prev) < tol:
            # 2. Distance and Azimuth Computation upon Convergence
            u2 = cos2_alpha * (a**2 - b**2) / (b**2)
            A = 1 + (u2 / 16384) * (4096 + u2 * (-768 + u2 * (320 - 175 * u2)))
            B = (u2 / 1024) * (256 + u2 * (-128 + u2 * (74 - 47 * u2)))
            
            delta_sigma = B * sin_sigma * (
                cos2_sigma_m + (B / 4) * (
                    cos_sigma * (-1 + 2 * cos2_sigma_m**2) -
                    (B / 6) * cos2_sigma_m * (-3 + 4 * sin_sigma**2) * (-3 + 4 * cos2_sigma_m**2)
                )
            )

            s = b * A * (sigma - delta_sigma)

            # Initial and final azimuths
            alpha1 = math.atan2(
                cosU2 * sin_lambda,
                cosU1 * sinU2 - sinU1 * cosU2 * cos_lambda
            )
            alpha2 = math.atan2(
                cosU1 * sin_lambda,
                -sinU1 * cosU2 + cosU1 * sinU2 * cos_lambda
            )

            # Normalize azimuths to [0, 360)
            az1 = (math.degrees(alpha1) + 360) % 360
            az2 = (math.degrees(alpha2) + 360) % 360

            return {
                'distance_m': s,
                'azimuth_initial_deg': az1,
                'azimuth_final_deg': az2,
                'converged': True
            }

    # If it fails to converge (e.g. nearly antipodal points)
    return {
        'distance_m': float('nan'),
        'azimuth_initial_deg': float('nan'),
        'azimuth_final_deg': float('nan'),
        'converged': False
    }

def get_altitude_km(rx_eci,ry_eci,rz_eci):
    """
    Get altitude from ECI point

    Parameters
    ----------
    rx_eci : float
        ECI x-coordinate
    ry_eci : float
        ECI y-coordinate
    rz_eci : float
        ECI z-coordinate

    Returns
    -------
    float
        altitude in kilometers
    """
    R = cct.eci_to_R(rx_eci,ry_eci,rz_eci)
    r = np.sqrt(rx_eci**2 + ry_eci**2 + rz_eci**2)
    return (r-R)/1000

def get_altitudes_km(rx_eci,ry_eci,rz_eci):
    """
    Get altitudes from ECI points

    Parameters
    ----------
    rx_eci : np.ndarray
        Vector of ECI x points
    ry_eci : TYPE
        Vector of ECI y points
    rz_eci : TYPE
        Vector of ECI z points

    Returns
    -------
    np.ndarray
        altitudes in kilometers
    """
    altitudes_km = np.empty(len(rx_eci))
    for idx in range(len(rx_eci)):
        altitude_km = get_altitude_km(rx_eci[idx],ry_eci[idx],rz_eci[idx])
        altitudes_km[idx] = altitude_km
    return altitudes_km

def get_ground_range_km(initial_lat_deg,initial_lon_deg,initial_ut1,current_r_eci,talo_sec):
    """
    Get the ground range from launch point given vehicle ECI position and corresponding time after launch

    Parameters
    ----------
    initial_lat_deg : float
        the launch point latitude
    initial_lon_deg : float
        the launch point longitude
    current_r_eci : np.ndarray
        Position vector in ECI [x, y, z] in meters.
    initial_ut1 : string
        the liftoff UT1 date-time string
    talo_sec : float
        time after liftoff in seconds

    Returns
    -------
    float
        ground range from launch point in kilometers
    """
    current_gmst = cct.get_current_gmst(initial_ut1,talo_sec)
    current_lat_deg,current_lon_deg,_ = cct.eci_to_lla(current_r_eci,current_gmst)
    vincenty_inv_dict = vincenty_inverse(initial_lat_deg,initial_lon_deg,current_lat_deg,current_lon_deg)
    return vincenty_inv_dict['distance_m']/1000

def get_ground_ranges_km(initial_lat_deg,initial_lon_deg,initial_ut1,rx_eci,ry_eci,rz_eci,talos_sec):
    """
    Get ground ranges from launch point given vehicle ECI positions and corresponding times after launch

    Parameters
    ----------
    initial_lat_deg : float
        the launch point latitude
    initial_lon_deg : float
        the launch point longitude
    rx_eci : np.ndarray
        Vector of ECI x points
    ry_eci : np.ndarray
        Vector of ECI y points
    rz_eci : np.ndarray
        Vector of ECI z points
    initial_ut1 : string
        the liftoff UT1 date-time string
    talos_sec : np.ndarray
        Vector of times since liftoff

    Returns
    -------
    np.ndarray
        ground ranges from launch point in meters
    """
    ground_ranges_km = np.empty(len(talos_sec))
    for idx,talo_sec in enumerate(talos_sec):
        current_gmst = cct.get_current_gmst(initial_ut1,talo_sec)
        current_r_eci = np.array([rx_eci[idx],ry_eci[idx],rz_eci[idx]])
        current_lat_deg,current_lon_deg,_ = cct.eci_to_lla(current_r_eci,current_gmst)
        vincenty_inv_dict = vincenty_inverse(initial_lat_deg,initial_lon_deg,current_lat_deg,current_lon_deg)
        ground_ranges_km[idx] = vincenty_inv_dict['distance_m']/1000
    return ground_ranges_km

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

    #s = get_ground_range(0,0,10,0)