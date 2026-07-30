# -*- coding: utf-8 -*-
"""
coordinate_conversion_tools.py

"""

import numpy as np
import geophysical_tools as gt

def ut1_to_jd(ut1_date_time):
    """
    Calculates the Julian Date given a UT1 date-time.
    The algorithm is based on the Gregorian calendar, which is valid for dates after 1582.
    Note: a Julian Day is a continuous count of days, while Julian Date is the exact day,
    which includes a fractional part.
    
    UT1 date-time format:
    Y-M-D H:Min:S
    """
    ut1_date, ut1_time = ut1_date_time.split(' ')
    ut1_date = ut1_date.split('-')
    Y = int(ut1_date[0])
    M = int(ut1_date[1])
    D = int(ut1_date[2])
    ut1_time = ut1_time.split(':')
    H = int(ut1_time[0])
    Min = int(ut1_time[1])
    S = int(ut1_time[2])
    
    if M <= 2:
        Y = Y - 1
        M = M + 12
    A = int(Y / 100)
    B = int(A / 4)
    C = 2 - A + B
    E = int(365.25 * (Y + 4716))
    F = int(30.6001 * (M + 1))
    day_fraction = (H + Min/60.0 + S/3600.0) / 24.0
    JD = C + D + day_fraction + E + F - 1524.5
    return JD

def jd_to_gmst(jd):
    """
    Calculates the Greenwich Mean Sidereal Time in radians
    Reference: "Astronomical Algorithms" by Jean Meeus.
    """
    #Calculate number of Julian centuries that have ellapsed since J2000.0 (2000-01-01 12:00:00 UT)
    jd_j2000 = 2451545.0
    T = (jd - jd_j2000)/36525
    
    #Calculate the sidereal angle in degrees using the IAU polynomial expression
    theta = 280.46061837 + 360.98564736629 * (jd - jd_j2000) + 0.000387933 * (T ** 2) - (T ** 3) / 38710000.0
    
    #Normalize angle to value in range [0,360)
    gmst_deg = theta % 360.0
    if gmst_deg < 0:
        gmst_deg += 360.0
        
    #Convert from degrees to radians
    gmst_rad = np.radians(gmst_deg)
    
    return gmst_rad

def jd_to_gmst_vec(jd_vec):
    """
    Calculates the Greenwich Mean Sidereal Times in radians for each Julian Date.
    Reference: "Astronomical Algorithms" by Jean Meeus.
    """
    #Calculate number of Julian centuries that have ellapsed since J2000.0 (2000-01-01 12:00:00 UT)
    jd_j2000_vec = np.ones(len(jd_vec)) * 2451545.0
    T_vec = (jd_vec - jd_j2000_vec)/36525
    
    #Calculate the sidereal angle in degrees using the IAU polynomial expression
    theta_vec = np.ones(len(jd_vec)) * 280.46061837 + np.ones(len(jd_vec)) * 360.98564736629 * (jd_vec - jd_j2000_vec) + 0.000387933 * (T_vec ** 2) - (T_vec ** 3) / 38710000.0
    
    #Normalize angle to value in range [0,360)
    gmst_deg_vec = theta_vec % 360.0
    for i in range(len(gmst_deg_vec)):
        if gmst_deg_vec[i] < 0:
            gmst_deg_vec[i] += 360.0
        
    #Convert from degrees to radians
    gmst_rad_vec = np.radians(gmst_deg_vec)
    
    return gmst_rad_vec

def ut1_to_gmst(ut1_date_time):
    '''
    Convert UT1 date-time to GMST in radians.
    UT1 date-time string format:
    Y-M-D H:Min:S

    Parameters
    ----------
    ut1_date_time : string
        the UT1 date-time

    Returns
    -------
    float
        GMST in radians
    '''
    jd = ut1_to_jd(ut1_date_time)
    gmst_rad = jd_to_gmst(jd)
    return gmst_rad

def get_current_gmst(liftoff_ut1,talo_sec):
    """
    Get the current Greenwich Mean Sidereal Time

    Parameters
    ----------
    liftoff_ut1 : string
        the UT1 date-time at liftoff
    talo_sec : float
        the time after lift off, in seconds

    Returns
    -------
    float
        the GMST in radians
    """
    liftoff_jd = ut1_to_jd(liftoff_ut1)
    talo_day = talo_sec/60/60/24
    current_jd = liftoff_jd + talo_day
    return jd_to_gmst(current_jd)

def get_gmsts(liftoff_ut1,talo_sec_array):
    """
    Get the Greenwich Mean Sidereal Times corresponding to the times in the 
    passed-in time after lift off array

    Parameters
    ----------
    liftoff_ut1 : string
        the UT1 date-time at liftoff
    talo_sec_array : np.array
        the array of times after lift off, in seconds

    Returns
    -------
    float
        the GMSTs in radians
    """
    liftoff_jd = ut1_to_jd(liftoff_ut1)
    liftoff_jd_array = np.ones(len(talo_sec_array)) * liftoff_jd
    talo_day_array = talo_sec_array/60/60/24
    current_jd_array = liftoff_jd_array + talo_day_array
    return jd_to_gmst_vec(current_jd_array)

def eci_to_ecef(r_eci: np.ndarray, v_eci: np.ndarray, theta: float):
    """
    Converts position and velocity from ECI to ECEF coordinates.
    Note that this doesn't account for the precession, nutation, or polar
    motion of Earth.

    Parameters:
    -----------
    r_eci : np.ndarray
        Position vector in ECI [x, y, z] in meters.
    v_eci : np.ndarray
        Velocity vector in ECI [vx, vy, vz] in meters/second.
    theta : float
        Greenwich Sidereal Time / Earth Rotation Angle in radians.

    Returns:
    --------
    r_ecef : np.ndarray
        Position vector in ECEF [x, y, z] in meters.
    v_ecef : np.ndarray
        Velocity vector in ECEF [vx, vy, vz] in meters/second.
    """
    cos_t = np.cos(theta)
    sin_t = np.sin(theta)

    # Rotation matrix around Z-axis by +theta
    R_z = np.array([[cos_t, sin_t, 0.0], [-sin_t, cos_t, 0.0], [0.0, 0.0, 1.0]])

    # 1. Transform Position
    r_ecef = R_z @ r_eci

    # 2. Transform Velocity
    # v_ecef = R_z * (v_eci - omega x r_eci)
    omega_vec = np.array([0.0, 0.0, gt.OMEGA_EARTH])
    v_transport = np.cross(omega_vec, r_eci)  # Transport velocity induced by Earth's rotation

    v_ecef = R_z @ (v_eci - v_transport)

    return r_ecef, v_ecef

def eci_to_ecef_vec(r_eci: np.ndarray, v_eci: np.ndarray, thetas: np.ndarray):
    """
    Converts position and velocity from ECI to ECEF coordinates.
    Note that this doesn't account for the precession, nutation, or polar
    motion of Earth.

    Parameters:
    -----------
    r_eci : np.ndarray
        An (N,3) NumPy array of position vectors in ECI [[x, y, z],...] in meters.
    v_eci : np.ndarray
        An (N,3) NumPy array of velocity vectors in ECI [vx, vy, vz] in meters/second.
    thetas : np.ndarray
        An (N,1) NumPy array of Greenwich Sidereal Times / Earth Rotation Angles in radians.

    Returns:
    --------
    r_ecef : np.ndarray
        An (N,3) NumPy array of position vectors in ECEF [x, y, z] in meters.
    v_ecef : np.ndarray
        An (N,3) NumPy array of position vectors in ECEF [vx, vy, vz] in meters/second.
    """
    cos_t = np.cos(thetas)
    sin_t = np.sin(thetas)
    
    # Rotation matrix around Z-axis by +theta
    #R_z = np.array([[cos_t, sin_t, 0.0], [-sin_t, cos_t, 0.0], [0.0, 0.0, 1.0]])

    # Transform Position
    #r_ecef = R_z @ r_eci
    
    # Directly apply the Z-axis rotation matrix row-by-row:
    # x_ecef =  cos(theta) * x + sin(theta) * y
    # y_ecef = -sin(theta) * x + cos(theta) * y
    # z_ecef =  z
    x_ecef = cos_t * r_eci[:, 0] + sin_t * r_eci[:, 1]
    y_ecef = -sin_t * r_eci[:, 0] + cos_t * r_eci[:, 1]
    z_ecef = r_eci[:, 2]

    r_ecef = np.column_stack((x_ecef, y_ecef, z_ecef))

    # Transform Velocity
    # v_ecef = R_z * (v_eci - omega x r_eci)
    # omega_vec = np.array([0.0, 0.0, gt.OMEGA_EARTH])
    # v_transport = np.cross(omega_vec, r_eci)  # Transport velocity induced by Earth's rotation
    #v_ecef = (R_z @ (v_eci - v_transport).T).T
    
    # Compute relative velocity in ECI (v_rel = v_eci - omega x r_eci)
    # Since omega = [0, 0, omega_e], omega x r = [-omega_e * y, omega_e * x, 0]
    v_rel = v_eci.copy()
    v_rel[:, 0] -= -gt.OMEGA_EARTH * r_eci[:, 1]
    v_rel[:, 1] -= gt.OMEGA_EARTH * r_eci[:, 0]

    # Apply the same Z-axis rotation to v_rel
    vx_ecef = cos_t * v_rel[:, 0] + sin_t * v_rel[:, 1]
    vy_ecef = -sin_t * v_rel[:, 0] + cos_t * v_rel[:, 1]
    vz_ecef = v_rel[:, 2]

    v_ecef = np.column_stack((vx_ecef, vy_ecef, vz_ecef))

    return r_ecef, v_ecef

def ecef_to_eci(r_ecef: np.ndarray, v_ecef: np.ndarray, theta: float):
    """
    Converts position and velocity from ECEF to ECI coordinates.
    Note that this doesn't account for the precession, nutation, or polar
    motion of Earth.

    Parameters:
    -----------
    r_ecef : np.ndarray
        Position vector in ECEF [x, y, z] in meters.
    v_ecef : np.ndarray
        Velocity vector in ECEF [vx, vy, vz] in meters/second.
    theta : float
        Greenwich Sidereal Time / Earth Rotation Angle in radians.

    Returns:
    --------
    r_eci : np.ndarray
        Position vector in ECI [x, y, z] in meters.
    v_eci : np.ndarray
        Velocity vector in ECI [vx, vy, vz] in meters/second.
    """
    cos_t = np.cos(theta)
    sin_t = np.sin(theta)

    # Inverse rotation matrix (transpose of R_z)
    R_z_inv = np.array(
        [[cos_t, -sin_t, 0.0], [sin_t, cos_t, 0.0], [0.0, 0.0, 1.0]]
    )

    # 1. Transform Position
    r_eci = R_z_inv @ r_ecef

    # 2. Transform Velocity
    # v_eci = R_z_inv * v_ecef + omega x r_eci
    omega_vec = np.array([0.0, 0.0, gt.OMEGA_EARTH])
    v_rel_eci = R_z_inv @ v_ecef
    v_transport = np.cross(omega_vec, r_eci)

    v_eci = v_rel_eci + v_transport

    return r_eci, v_eci

def ecef_to_eci_vec(r_ecef: np.ndarray, v_ecef: np.ndarray, thetas: float):
    """
    Converts position and velocity from ECEF to ECI coordinates.
    Note that this doesn't account for the precession, nutation, or polar
    motion of Earth.

    Parameters:
    -----------
    r_ecef : np.ndarray
        An (N,3) NumPy array of position vectors in ECEF [[x, y, z],...] in meters.
    v_ecef : np.ndarray
        An (N,3) NumPy array of velocity vectors in ECEF [vx, vy, vz] in meters/second.
    theta : np.ndarray
        Greenwich Sidereal Time / Earth Rotation Angle in radians.

    Returns:
    --------
    r_eci : np.ndarray
        An (N,3) NumPy array of position vectors in ECI [[x, y, z],...] in meters.
    v_eci : np.ndarray
        An (N,3) NumPy array of position vectors in ECI [[vx, vy, vz],...] in meters/second.
    """

    cos_t = np.cos(thetas)
    sin_t = np.sin(thetas)

    # Inverse rotation matrix (transpose of R_z)
    #R_z_inv = np.array([[cos_t, -sin_t, 0.0], [sin_t, cos_t, 0.0], [0.0, 0.0, 1.0]])

    # Transform Position
    #multiplying 3x3 rotation matrix with 3xN matrix to get 3xN matrix, which is transposed to get Nx3 matrix
    #r_eci = (R_z_inv @ r_ecef.T).T 
    
    x_eci = cos_t * r_ecef[:,0] - sin_t * r_ecef[:,1]
    y_eci = sin_t * r_ecef[:,0] + cos_t * r_ecef[:,1]
    z_eci = r_ecef[:,2]

    r_eci = np.column_stack((x_eci,y_eci,z_eci))

    # Transform Velocity
    # v_eci = R_z_inv * v_ecef + omega x r_eci
    #omega_vec = np.array([0.0, 0.0, gt.OMEGA_EARTH])
    #v_rel_eci = (R_z_inv @ v_ecef.T).T
    #v_transport = np.cross(omega_vec, r_eci) #this compute omega_vec x r_eci for every row in r_eci

    vx_rel_eci = cos_t * v_ecef[:,0] - sin_t * v_ecef[:,1]
    vy_rel_eci = sin_t * v_ecef[:,0] + cos_t * v_ecef[:,1]
    vz_rel_eci = v_ecef[:,2]
    
    vx_transport = -gt.OMEGA_EARTH * r_eci[:,1]
    vy_transport = gt.OMEGA_EARTH * r_eci[:,0]

    v_eci = np.column_stack((vx_rel_eci+vx_transport,vy_rel_eci+vy_transport,vz_rel_eci))

    return r_eci, v_eci

def eci_to_lla(r_eci,theta, degrees: bool = True) -> tuple[float, float, float]:
    """
    Converts a single ECI point to Latitude, Longitude, and Altitude

    Parameters
    ----------
    r_eci : np.ndarray
        Position vector in ECI [x, y, z] in meters.
    theta : float
        the Greenwich Mean Sidereal Time in radians

    Returns
    -------
    lat : float
        the latitude
    lon : float
        the longitude
    alt : float
        the altitude, in meters
    """
    
    r_ecef,_ = eci_to_ecef(r_eci, [0,0,0], theta)
    x = r_ecef[0]
    y = r_ecef[1]
    z = r_ecef[2]
    lat,lon,alt = ecef_to_lla(x,y,z,degrees=False)
    if degrees:
        lat = np.degrees(lat)
        lon = np.degrees(lon)
    return lat,lon,alt

def eci_to_lla_vec(r_eci: np.ndarray,thetas: np.ndarray, degrees: bool = True) -> tuple[float, float, float]:
    """
    Converts ECI points to Latitude, Longitude, and Altitude points

    Parameters
    ----------
    r_eci : np.ndarray
        An (N,3) NumPy array of position vectors in ECI [[x, y, z],...] in meters.
    thetas : float
        the Greenwich Mean Sidereal Times in radians

    Returns
    -------
    lat : np.ndarray
        An (N,3) NumPy array of latitudes
    lon : np.ndarray
        An (N,3) NumPy array of longitudes
    alt : np.ndarray
        An (N,3) NumPy array of altitudes, in meters
    """
    
    r_ecef,_ = eci_to_ecef_vec(r_eci, np.zeros((len(r_eci), 3), dtype=np.float64), thetas)
    lla = ecef_to_lla_vec(r_ecef,degrees=False)
    if degrees==True:
        lla[:,0] = np.degrees(lla[:,0])
        lla[:,1] = np.degrees(lla[:,1])
    return lla

def lla_to_ecef(lat_deg: float, lon_deg: float, alt_m: float) -> tuple[float, float, float]:
    """
    Converts Geodetic coordinates (LLA) to Earth-Centered, Earth-Fixed (ECEF) Cartesian coordinates.

    Parameters
    ----------
    lat_deg (float): 
        Latitude in degrees (-90 to 90)
    lon_deg (float): 
        Longitude in degrees (-180 to 180)
    alt_m (float): 
        Altitude above ellipsoid in meters

    Returns
    -------
        Tuple[float, float, float]: ECEF coordinates (X, Y, Z) in meters
    """
    # Convert latitude and longitude from degrees to radians
    phi = np.radians(lat_deg)
    lam = np.radians(lon_deg)

    # Prime vertical radius of curvature
    sin_phi = np.sin(phi)
    cos_phi = np.cos(phi)
    sin_lam = np.sin(lam)
    cos_lam = np.cos(lam)

    # Calculate the radius of curvature in the prime vertical
    N = gt.R_eq / np.sqrt(1.0 - gt.e2 * (sin_phi ** 2))

    # Calculate Cartesian coordinates
    x = (N + alt_m) * cos_phi * cos_lam
    y = (N + alt_m) * cos_phi * sin_lam
    z = (N * (1.0 - gt.e2) + alt_m) * sin_phi

    return x, y, z

def lla_to_ecef_vec(lla_pts: np.ndarray):
    """
    Converts Geodetic coordinates (LLA) to Earth-Centered, Earth-Fixed (ECEF) Cartesian coordinates.
    
    The input is a (N, 3) NumPy array of LLA coordinates [[lat, lon, alt], ...].
    Returns (N, 3) array of [[x, y, z], ...].
    """
    lat_deg, lon_deg, alt_m = lla_pts[:, 0], lla_pts[:, 1], lla_pts[:, 2]
    
    phi = np.radians(lat_deg)
    lam = np.radians(lon_deg)

    sin_phi = np.sin(phi)
    cos_phi = np.cos(phi)

    N = gt.R_eq / np.sqrt(1.0 - gt.e2 * (sin_phi ** 2))

    x = (N + alt_m) * cos_phi * np.cos(lam)
    y = (N + alt_m) * cos_phi * np.sin(lam)
    z = (N * (1.0 - gt.e2) + alt_m) * sin_phi

    return np.column_stack((x, y, z))

def lla_to_eci(lat_deg: float,lon_deg: float,alt_m: float,v_rel: np.ndarray,theta: float):
    """
    Converts a single LLA point to an ECI Cartesian position and velocity vector

    Parameters
    ----------
    lat_deg : float
        the latitude in degrees
    lon_deg : float
        the longitude in degrees
    alt_m : float
        the altitude in meters
    v_rel : a (1, 3) np.ndarray
        the velocity vector [vx,vy,vz]
    theta : float
        the GMST in radians

    Returns
    -------
    r_eci : a (1, 3) np.ndarray
        the position vector [x,y,z]
    v_eci : a (1, 3) np.ndarray
        the velocity vector [vx,vy,vz]
    """
    
    x_ecef,y_ecef,z_ecef = lla_to_ecef(lat_deg,lon_deg,alt_m)
    r_ecef = [x_ecef,y_ecef,z_ecef]
    r_eci,v_eci = ecef_to_eci(r_ecef, v_rel, theta)
    return r_eci,v_eci

def lla_to_eci_vec(lla_pts:np.ndarray,v_ecef: np.ndarray,thetas: float):
    """
    Converts a single LLA point to an ECI Cartesian position and velocity vector

    Parameters
    ----------
    lla_pts: an (N, 3) np.ndarray
        the LLA points [[lat,lon,alt],...]
    v_ecef : an (N, 3) np.ndarray
        the relative velocity vectors [[vx,vy,vz],...] measured in ECEF frame
    thetas : an (N,1) np.ndarray
        the GMSTs in radians

    Returns
    -------
    r_eci : a (N, 3) np.ndarray
        the position vector [x,y,z]
    v_eci : a (N, 3) np.ndarray
        the velocity vector [vx,vy,vz]
    """
    
    r_ecef = lla_to_ecef_vec(lla_pts)
    r_eci,v_eci = ecef_to_eci_vec(r_ecef, v_ecef, thetas)
    return r_eci,v_eci

def ecef_to_lla(x: float, y: float, z: float, degrees: bool = True) -> tuple[float, float, float]:
    """
    Converts a single ECEF Cartesian coordinate (X, Y, Z in meters) to 
    a geodetic LLA coordinate (Latitude, Longitude, Altitude).
    
    Uses Bowring's closed-form algorithm (1985).
    """
    # Step 1: Distance from Z-axis
    p = np.hypot(x, y)
    
    # Handle origin/center of Earth edge case
    if p == 0 and z == 0:
        return 0.0, 0.0, -gt.R_eq

    # Step 2: Longitude
    lon = np.atan2(y, x)
    
    # Step 3: Auxiliary parametric angle theta
    theta = np.atan2(gt.R_eq * z, gt.R_p * p)
    
    # Step 4: Latitude
    sin_theta = np.sin(theta)
    cos_theta = np.cos(theta)
    
    num = z + gt.eprime2 * gt.R_p * (sin_theta ** 3)
    den = p - gt.e2 * gt.R_eq * (cos_theta ** 3)
    lat = np.atan2(num, den)
    
    # Step 5: Altitude
    sin_lat = np.sin(lat)
    cos_lat = np.cos(lat)
    
    # Radius of curvature in the prime vertical
    N = gt.R_eq / np.sqrt(1.0 - gt.e2 * (sin_lat ** 2))
    
    # Pole boundary handling to prevent divide-by-zero near cos(lat) == 0
    if abs(cos_lat) > 1e-6:
        alt = (p / cos_lat) - N
    else:
        alt = (abs(z) / sin_lat) - N * (1.0 - gt.e2)
        
    if degrees:
        return np.degrees(lat), np.degrees(lon), alt
    
    return lat, lon, alt

def ecef_to_lla_vec(ecef_pts: np.ndarray, degrees: bool = True) -> np.ndarray:
    """
    Converts ECEF Cartesian coordinates (X, Y, Z in meters) to 
    Geodetic LLA coordinates (Latitude, Longitude, Altitude).
    
    Uses Bowring's closed-form algorithm (1985).
    
    Vectorized conversion for (N, 3) NumPy array of ECEF coordinates [[X, Y, Z], ...].
    Returns (N, 3) array of [[Lat, Lon, Alt], ...].
    """
    x, y, z = ecef_pts[:, 0], ecef_pts[:, 1], ecef_pts[:, 2]
    
    p = np.hypot(x, y)
    lon = np.arctan2(y, x)
    theta = np.arctan2(gt.R_eq * z, gt.R_p * p)
    
    sin_theta = np.sin(theta)
    cos_theta = np.cos(theta)
    
    num = z + gt.eprime2 * gt.R_p * (sin_theta ** 3)
    den = p - gt.e2 * gt.R_eq * (cos_theta ** 3)
    lat = np.arctan2(num, den)
    
    sin_lat = np.sin(lat)
    cos_lat = np.cos(lat)
    
    N = gt.R_eq / np.sqrt(1.0 - gt.e2 * (sin_lat ** 2))
    
    # Use general formula, mask poles if necessary
    alt = np.where(np.abs(cos_lat) > 1e-6,
                   (p / cos_lat) - N,
                   (np.abs(z) / np.abs(sin_lat)) - N * (1.0 - gt.e2))
    
    if degrees:
        lat = np.degrees(lat)
        lon = np.degrees(lon)
        
    return np.column_stack((lat, lon, alt))

def lat_lon_to_R(lat,lon):
    '''
    Returns radius of WGS84 Earth at the given latitude and longitude.
    
    Parameters
    ----------
    lat : float
        latitude in radians
    lon : float
        longitude in radians

    Returns
    -------
    float
        radius of WGS84 Earth
    '''
    return np.sqrt(((gt.R_eq**2*np.cos(lat))**2 + (gt.R_p**2 * np.sin(lat))**2)/((gt.R_eq*np.cos(lat))**2 + (gt.R_p * np.sin(lat))**2))

def eci_to_R(x_eci,y_eci,z_eci):
    '''
    Get radius of Earth directly beneath the ECI point. This assumes that the
    coordinate transformation from ECI to ECEF is a rotation about the z-axis only
    (i.e. precession, nutation, and polar motion are not modeled in the rotation matrix)
    Given this assumption, the ECI point can be directly used as the ECEF point in the 
    conversion from ECEF to LLA, since we need latitude and it doesn't change with
    rotation about the z-axis.
    
    Parameters
    ----------
    x_eci : float
        the ECI x-coordinate
    y_eci : float
        the ECI y-coordinate
    z_eci : float
        the ECI z-coordinate

    Returns
    -------
    float
        radius of Earth directly beneath the ECI point, in meters
    '''
    lat,_,_ = ecef_to_lla(x_eci,y_eci,z_eci,degrees=False)
    return np.sqrt(((gt.R_eq**2*np.cos(lat))**2 + (gt.R_p**2 * np.sin(lat))**2)/((gt.R_eq*np.cos(lat))**2 + (gt.R_p * np.sin(lat))**2))

if __name__ == "__main__":
    
    r_eci_init = np.array([5000000.0, 3000000.0, 4000000.0])  # meters
    v_eci_init = np.array([-3000.0, 6000.0, 1000.0])  # m/s
    theta = np.radians(45)
    
    r_ecef1,v_ecef1 = eci_to_ecef(r_eci_init, v_eci_init, theta)
    r_ecef_vec1 = np.array([r_ecef1,r_ecef1])
    v_ecef_vec1 = np.array([v_ecef1,v_ecef1])
    
    r_eci1, v_eci1 = ecef_to_eci(r_ecef1, v_ecef1, theta)
    
    r_eci_vec1,v_eci_vec1 = ecef_to_eci_vec(r_ecef_vec1, v_ecef_vec1, theta)
    
    '''
    # Example state in ECI (Low Earth Orbit satellite)
    r_eci_init = np.array([5000000.0, 3000000.0, 4000000.0])  # meters
    v_eci_init = np.array([-3000.0, 6000.0, 1000.0])  # m/s

    # Example Earth Rotation Angle (theta) = 45 degrees in radians
    theta_gst = np.radians(45.0)

    # --- Convert ECI -> ECEF ---
    r_ecef, v_ecef = eci_to_ecef(r_eci_init, v_eci_init, theta_gst)

    print("=== ECI to ECEF Conversion ===")
    print(f"ECEF Position (m) : {r_ecef}")
    print(f"ECEF Velocity (m/s): {v_ecef}\n")

    # --- Convert ECEF -> ECI (Round-Trip Test) ---
    r_eci_reconstructed, v_eci_reconstructed = ecef_to_eci(
        r_ecef, v_ecef, theta_gst
    )

    print("=== Round-Trip Validation (ECEF -> ECI) ===")
    print(f"Original ECI Position : {r_eci_init}")
    print(f"Reconstructed ECI Pos : {r_eci_reconstructed}")
    print(f"Position Error (m)    : {np.linalg.norm(r_eci_init - r_eci_reconstructed):.2e}\n")

    print(f"Original ECI Velocity : {v_eci_init}")
    print(f"Reconstructed ECI Vel : {v_eci_reconstructed}")
    print(f"Velocity Error (m/s)  : {np.linalg.norm(v_eci_init - v_eci_reconstructed):.2e}")


    # Test for 2026-07-23 19:44:00
    ut1_date_time1 = '2026-07-23 19:44:00'
    jd1 = ut1_to_jd(ut1_date_time1)
    print(f'The Julian Date of {ut1_date_time1} is {jd1}')
    
    # Test for 2000-01-01 12:00:00
    ut1_date_time2 = '2000-01-01 12:00:00'
    jd2 = ut1_to_jd(ut1_date_time2)
    print(f'The Julian Date of {ut1_date_time2} is {jd2}')
    gmst2 = jd_to_gmst(jd2)
    print(f'The GMST at {ut1_date_time2} is {round(np.degrees(gmst2))} degrees')
    
    # Test for 2000-03-20 12:00:00
    ut1_date_time3 = '2026-07-24 12:00:00'
    jd3 = ut1_to_jd(ut1_date_time3)
    print(f'The Julian Date of {ut1_date_time3} is {jd3}')
    gmst3 = jd_to_gmst(jd3)
    print(f'The GMST at {ut1_date_time3} is {round(np.degrees(gmst3),4)} degrees')
    '''
    '''
    # Test Point: Space Needle area (approx 47.6205° N, -122.3493° W, 100m alt)
    # Corresponding ECEF coordinates:
    x_test, y_test, z_test = -2280521.14, -3632738.92, 4689408.26
    
    lat, lon, alt = ecef_to_lla(x_test, y_test, z_test)
    
    print("--- Single Point Result ---")
    print(f"ECEF: ({x_test}, {y_test}, {z_test})")
    print(f"Latitude:  {lat:.6f}°")
    print(f"Longitude: {lon:.6f}°")
    print(f"Altitude:  {alt:.2f} m")
    
    # Batch Test
    points = np.array([
        [x_test, y_test, z_test],
        [6378137.0, 0.0, 0.0],  # Prime Meridian on Equator
        [0.0, 0.0, gt.R_p]     # North Pole surface
    ])
    
    lla_batch = ecef_to_lla_vec(points)
    print("\n--- Vectorized Batch Results (Lat, Lon, Alt) ---")
    print(np.round(lla_batch, 4))


    # Example: Kennedy Space Center Launch Pad 39A
    lat = 28.608389   # degrees North
    lon = -80.604333  # degrees East
    alt = 0.0         # meters

    x, y, z = lla_to_ecef(lat, lon, alt)
    print(f"Input LLA  : Lat = {lat}°, Lon = {lon}°, Alt = {alt} m")
    print(f"Output ECEF: X = {x:.3f} m, Y = {y:.3f} m, Z = {z:.3f} m")
    
    ecef_batch = lla_to_ecef_vec(lla_batch)
'''