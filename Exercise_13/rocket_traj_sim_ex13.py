# -*- coding: utf-8 -*-
"""
rocket_traj_sim_ex13.py

Description:
    The current rocket thrust equation uses mass flow rate, m_dot, and exhaust velocity, v_exhaust:
    T = m_dot * v_exhaust
    
    Update the rocket thrust equation to take into account the change in atmospheric pressure, that is, account for the pressure thrust:
    T = m_dot * v_exhaust + (p_exhaust - p(h)) * A_exit
    
    Refactor the atmospheric model to use the U.S. Standard Atmosphere 1976 specification, then calculate 
    pressure as a function  of height, p(h). Additionally, update the density, temperature, and speed of sound calculations to use this more
    accurate atmospheric model. Finally, update the RocketStage class to have the p_exhaust and A_exit properties required for computing thrust.
    
    References:
        - U.S. Standard Atmosphere model, 0km to 86km: 'Calculating Atmospheric Conditions (Temperature, Pressure, Air Density, and Speed of Sound)
          Using C++' by Robert Yager
        - U.S. Standard Atmosphere model, 86km to 1000km: http://www.braeunig.us/space/atmmodel.htm

Created by: Bradley Canty, 2026/07/21

TO DO: 

"""

import numpy as np
from Rocket import Rocket, RocketStage

v2_s1 = RocketStage(
    name='stage 1',
    m_empty=3960,
    m_propellant=8840,
    m_dot=129.4,
    v_exhaust=2050,
    p_exhaust = 101325, #Pa, engine pressure at nozzle exit plane
    A_exit = 0.43, #m^2, total engine nozzle exit area
    stage_diameter_m=1.65, #m^2
    coast_duration_sec = 0, #sec, coast duration (occurs before engine start)
    t_steer_start_rel = 3, #sec, steering start time relative to stage's engine start, If first stage, this is the pitch-over time
    steer_duration_sec = 10, #sec, duration to apply thrust offset by angle eta from the velocity vector
    target_fpa_at_bo_deg = 45, #deg, target flight path angle at stage burnout
    )

v2 = Rocket(
    {
        'name'               : 'V2',
		'date0'              : '2020-01-01 00:00:00', #initial date-time
		'frame'              : 'J2000',      #standard for earth-centered propagation
		'ode_solver'         : 'rk4',
        'lla0'               : np.array([0,0,0]), #Latitude [degrees], Longitude [degrees], altitude [meters]
        'v0_rel'             : np.array([0,0,0]), #vx [m/s], vy [m/s], vz [m/s]
		'stages'             : [v2_s1],
        'tspan'              : 1000,
        'dt'                 : .5,
        'm_payload'          : 1000,
        'payload_separates'  : False,
        'launch_azimuth_deg' : 90, #since this is a 2D simulation, set this to either 90 or -90
	})

mm3_s1 = RocketStage(
    name='stage 1',
    m_empty=2060, #kg
    m_propellant=20950, #kg
    m_dot=392, #kg/s
    v_exhaust=2635, #m/s
    p_exhaust = 55000, #Pa, engine pressure at nozzle exit plane
    A_exit = 0.445,    #m^2, total engine nozzle exit area
    stage_diameter_m=1.68, #m
    coast_duration_sec = 0, #sec, coast duration (occurs before engine start)
    t_steer_start_rel = 3, #sec, steering start time relative to stage's engine start
    steer_duration_sec = 15, #sec, duration to apply thrust offset by angle eta from the velocity vector; set to -1 to steer for the engine burn time
    target_fpa_at_bo_deg = 45, #deg, target flight path angle at stage burnout
    )

mm3_s2 = RocketStage(
    name='stage 2',
    m_empty=680, #kg
    m_propellant=6430, #kg
    m_dot=94, #kg/s
    v_exhaust=2824, #m/s
    p_exhaust = 100000, #Pa, engine pressure at nozzle exit plane
    A_exit = 0.639,    #m^2, total engine nozzle exit area
    stage_diameter_m=1.32, #m
    coast_duration_sec = 0, #sec, coast duration (occurs before engine start)
    t_steer_start_rel = 0, #sec, steering start time relative to stage's engine start
    steer_duration_sec = -1, #sec, duration to apply thrust offset by angle eta from the velocity vector
    target_fpa_at_bo_deg = 28, #deg, target flight path angle at stage burnout
    )

mm3_s3 = RocketStage(
    name='stage 3',
    m_empty=350, #kg
    m_propellant=3310, #kg
    m_dot=51, #kg/s
    v_exhaust=2844, #m/s
    p_exhaust = 38000, #Pa, engine pressure at nozzle exit plane
    A_exit = 0.4,    #m^2, total engine nozzle exit area
    stage_diameter_m=1.32, #m
    coast_duration_sec = 0, #sec, coast duration (occurs before engine start)
    t_steer_start_rel = 0, #sec, steering start time relative to stage's engine start
    steer_duration_sec = -1, #sec, duration to apply thrust offset by angle eta from the velocity vector
    target_fpa_at_bo_deg = 20, #deg, target flight path angle at stage burnout
    )

mm3 = Rocket(
    {
        'name'               : 'Minuteman 3',
		'date0'              : '2020-01-01', #initial date
		'et0'                : None,         #initial ephemeris time, i.e., seconds since J2000
		'frame'              : 'J2000',      #standard for earth-centered propagation
		'ode_solver'         : 'rk4',
        'lla0'               : np.array([0,0,0]), #Latitude [degrees], Longitude [degrees], altitude [meters]
        'v0_rel'             : np.array([0,0,0]), #vx [m/s], vy [m/s], vz [m/s]
		'stages'             : [mm3_s1,mm3_s2,mm3_s3],
        'tspan'              : 2*60*60,
        'dt'                 : .5,
        'm_payload'          : 1250, #kg, 750 kg mass of PBV and shroud + 500 kg mass of 1 RV
        'payload_separates'  : True,
        'launch_azimuth_deg' : 90,
	})

f9_s1 = RocketStage(
    name='stage 1',
    m_empty=25600, #kg
    m_propellant=395700, #kg
    m_dot=2750, #kg/s
    v_exhaust=2800, #m/s
    p_exhaust = 100000, #Pa, engine pressure at nozzle exit plane
    A_exit = 6.24,    #m^2, total engine nozzle exit area
    stage_diameter_m=3.7, #m
    coast_duration_sec = 0, #sec, coast duration (occurs before engine start)
    t_steer_start_rel = 5, #sec, steering start time relative to stage's engine start
    steer_duration_sec = 0.5, #sec, duration to apply thrust offset by angle eta from the velocity vector; set to -1 to steer for the engine burn time
    target_fpa_at_bo_deg = 45, #deg, target flight path angle at stage burnout
    )

f9_s2 = RocketStage(
    name='stage 2',
    m_empty=3900, #kg
    m_propellant=92670, #kg
    m_dot=287, #kg/s
    v_exhaust=3410, #m/s
    p_exhaust = 1000, #Pa, engine pressure at nozzle exit plane
    A_exit = 5.9,    #m^2, total engine nozzle exit area
    stage_diameter_m=3.7, #m
    coast_duration_sec = 11, #sec, coast duration (occurs before engine start)
    t_steer_start_rel = 0, #sec, steering start time relative to stage's engine start
    steer_duration_sec = -1, #sec, duration to apply thrust offset by angle eta from the velocity vector
    target_fpa_at_bo_deg = 0, #deg, target flight path angle at stage burnout
    )

f9 = Rocket(
    {
        'name'               : 'Falcon 9',
		'date0'              : '2020-01-01', #initial date
		'et0'                : None,         #initial ephemeris time, i.e., seconds since J2000
		'frame'              : 'J2000',      #standard for earth-centered propagation
		'ode_solver'         : 'rk4',
        'lla0'               : np.array([0,0,0]), #Latitude [degrees], Longitude [degrees], altitude [meters]
        'v0_rel'             : np.array([0,0,0]), #vx [m/s], vy [m/s], vz [m/s]
		'stages'             : [f9_s1,f9_s2],
        'tspan'              : 2*60*60,
        'dt'                 : .5,
        'm_payload'          : 17000, #kg
        'payload_separates'  : True,
        'launch_azimuth_deg' : 90,
	})

v2.generate_all_plots(output_img=True,output_gif=True,gif_frame_count=50)
mm3.generate_all_plots(output_img=True,output_gif=True,gif_frame_count=150)
f9.generate_all_plots(output_img=True,output_gif=True,gif_frame_count=150)

#v2.plot_2d_traj(rel_to_launch_point=False,output_img=True,output_gif=False,gif_frame_count=150)
#v2.plot_mass_vs_time(output_img=True,show_flight_phases=True)
#mm3.plot_mass_vs_time(output_img=True,show_flight_phases=True)
#v2.plot_mass_vs_time(output_img=True,show_flight_phases=True)






