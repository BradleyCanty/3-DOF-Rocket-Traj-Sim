# -*- coding: utf-8 -*-
"""
rocket_traj_sim_ex15.py

Description:
    Extend the simulation from 2 dimensions to 3 dimensions, and plot a
    rocket trajectory on a 3D globe.

Created by: Bradley Canty, 2026/07/28
"""

import numpy as np
from Rocket import Rocket, RocketStage

v2_s1 = RocketStage(
    name='stage 1',
    m_empty=3960,
    m_propellant=8840,
    m_dot=130,
    v_exhaust=2050,
    p_exhaust = 83000, #Pa, exhaust pressure at nozzle exit plane
    A_exit = 0.44, #m^2, total engine nozzle exit area
    stage_diameter_m=1.65, #m^2
    coast_duration_sec = 0, #sec, coast duration (occurs before engine start)
    t_steer_start_rel = 3, #sec, steering start time relative to stage's engine start, If first stage, this is the pitch-over time
    steer_duration_sec = 10, #sec, duration to apply thrust offset by angle eta from the velocity vector
    target_fpa_at_bo_deg = 40, #deg, target flight path angle at stage burnout
    )

v2 = Rocket(
    {
        'name'               : 'V2',
		'ut1datetime0'       : '2020-04-01 00:00:00', #initial UT1 date-time; format: year-month-day hours:minutes:seconds
		'ode_solver'         : 'rk4',
        'lla0'               : np.array([45,0,0]), #Latitude [degrees], Longitude [degrees], altitude [meters]
        'v0_rel'             : np.array([0,0,0]), #vx [m/s], vy [m/s], vz [m/s]
		'stages'             : [v2_s1],
        'tspan'              : 1000,
        'dt'                 : .5,
        'm_payload'          : 1000,
        'payload_separates'  : False,
        'launch_azimuth_deg' : 90,
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
    target_fpa_at_bo_deg = 40, #deg, target flight path angle at stage burnout
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
		'ut1datetime0'       : '2020-01-01 00:00:00', #initial UT1 date-time; format: year-month-day hours:minutes:seconds
		'et0'                : None,         #initial ephemeris time, i.e., seconds since J2000
		'ode_solver'         : 'rk4',
        'lla0'               : np.array([45,0,0]), #Latitude [degrees], Longitude [degrees], altitude [meters]
        'v0_rel'             : np.array([0,0,0]), #vx [m/s], vy [m/s], vz [m/s]
		'stages'             : [mm3_s1,mm3_s2,mm3_s3],
        'tspan'              : 2*60*60,
        'dt'                 : .5,
        'm_payload'          : 2250, #kg, 750 kg mass of PBV and shroud + 1500 kg mass of 3 RVs (500 kg per RV)
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
		'ut1datetime0'       : '2020-01-01 00:00:00', #initial UT1 date-time; format: year-month-day hours:minutes:seconds
		'et0'                : None,         #initial ephemeris time, i.e., seconds since J2000
		'ode_solver'         : 'rk4',
        'lla0'               : np.array([28.57,-80.65,0]), #Latitude [degrees], Longitude [degrees], altitude [meters]
        'v0_rel'             : np.array([0,0,0]), #vx [m/s], vy [m/s], vz [m/s]
		'stages'             : [f9_s1,f9_s2],
        'tspan'              : 2*60*60,
        'dt'                 : .5,
        'm_payload'          : 17000, #kg
        'payload_separates'  : True,
        'launch_azimuth_deg' : 90,
	})

v2.generate_all_plots(output_img=True)
mm3.generate_all_plots(output_img=True)
f9.generate_all_plots(output_img=True)

f9.show_traj_on_globe()