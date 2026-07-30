# -*- coding: utf-8 -*-
"""
rocket_traj_sim_ex11.py

Description:
    Continuing from exercise 10, extend the simulation to have a rocket class which
    instantiates a rocket object for each rocket being simulated. Organize the plotting
    code such that plots are generated via method calls on rocket objects. Additionally, 
    extend the code to handle rockets with multiple stages.
    
    Finally, produce the same plots from exercise 9, that is, plot...
    1) the rocket trajectory (i.e. altitude versus downrange distance)
    2) vertical distance vs horizontal distance from launch site across limb of earth
    3) flight path angle vs time
    4) position, velocity, and acceleration vs time
    5) dynamic pressure vs time
    6) mass vs time
    7) Mach number and drag coefficient vs time
    8) a GIF of the trajectory
    
    but for...
    1) V2 rocket (single stage, suborbital)
    2) LGM-30B Minuteman I ICBM (three stage, suborbital)
    3) Falcon 9 (two stage, orbital).

Created by: Bradley Canty, 2026/07/07

TO DO: 

"""

import numpy as np
from Rocket import Rocket, RocketStage

v2_s1 = RocketStage(
    name='stage 1',
    m_empty=3000,
    m_propellant=8840,
    m_dot=129.4,
    v_exhaust=2050,
    stage_diameter_m=1.65,
    coast_duration_sec = 0, #sec, coast duration (occurs before engine start)
    t_steer_start_rel = 3, #sec, steering start time relative to stage's engine start, If first stage, this is the pitch-over time
    steer_duration_sec = 10, #sec, duration to apply thrust offset by angle eta from the velocity vector
    target_fpa_at_bo_deg = 45, #deg, target flight path angle at stage burnout
    )

v2 = Rocket(
    {
        'name'               : 'V2',
		'date0'              : '2020-01-01', #initial date
		'et0'                : None,         #initial ephemeris time, i.e., seconds since J2000
		'frame'              : 'J2000',      #standard for earth-centered propagation
		'ode_solver'         : 'rk4',
        'lla0'               : np.array([0,0,0]), #Latitude [degrees], Longitude [degrees], altitude [meters]
        'v0_rel'             : np.array([0,0,0]), #vx [m/s], vy [m/s], vz [m/s]
		'stages'             : [v2_s1],
        'tspan'              : 1000,
        'dt'                 : .5,
        'm_payload'          : 1000,
        'payload_separates'  : False,
        'launch_azimuth_deg' : 90,   #Normally changes the sign on eta, right now does nothing...
	})

mm3_s1 = RocketStage(
    name='stage 1',
    m_empty=2292, #kg
    m_propellant=20785, #kg
    m_dot=343, #kg/s
    v_exhaust=2635, #m/s
    stage_diameter_m=1.68, #m
    coast_duration_sec = 0, #sec, coast duration (occurs before engine start)
    t_steer_start_rel = 3, #sec, steering start time relative to stage's engine start
    steer_duration_sec = 15, #sec, duration to apply thrust offset by angle eta from the velocity vector; set to -1 to steer for the engine burn time
    target_fpa_at_bo_deg = 40, #deg, target flight path angle at stage burnout
    )

mm3_s2 = RocketStage(
    name='stage 2',
    m_empty=795, #kg
    m_propellant=6237, #kg
    m_dot=94.8, #kg/s
    v_exhaust=2824, #m/s
    stage_diameter_m=1.32, #m
    coast_duration_sec = 0, #sec, coast duration (occurs before engine start)
    t_steer_start_rel = 0, #sec, steering start time relative to stage's engine start
    steer_duration_sec = -1, #sec, duration to apply thrust offset by angle eta from the velocity vector
    target_fpa_at_bo_deg = 20, #deg, target flight path angle at stage burnout
    )

mm3_s3 = RocketStage(
    name='stage 3',
    m_empty=400, #kg
    m_propellant=3200, #kg
    m_dot=46.1, #kg/s
    v_exhaust=2844, #m/s
    stage_diameter_m=1.32, #m
    coast_duration_sec = 0, #sec, coast duration (occurs before engine start)
    t_steer_start_rel = 0, #sec, steering start time relative to stage's engine start
    steer_duration_sec = -1, #sec, duration to apply thrust offset by angle eta from the velocity vector
    target_fpa_at_bo_deg = 8, #deg, target flight path angle at stage burnout
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
        'm_payload'          : 1000, #kg
        'payload_separates'  : True,
        'launch_azimuth_deg' : 90,
	})

f9_s1 = RocketStage(
    name='stage 1',
    m_empty=22200, #kg
    m_propellant=411000, #kg, CHECK THIS
    m_dot=2750, #kg/s
    v_exhaust=2800, #m/s
    stage_diameter_m=3.7, #m
    coast_duration_sec = 0, #sec, coast duration (occurs before engine start)
    t_steer_start_rel = 5, #sec, steering start time relative to stage's engine start
    steer_duration_sec = 0.5, #sec, duration to apply thrust offset by angle eta from the velocity vector; set to -1 to steer for the engine burn time
    target_fpa_at_bo_deg = 45, #deg, target flight path angle at stage burnout
    )

f9_s2 = RocketStage(
    name='stage 2',
    m_empty=4500, #kg
    m_propellant=111500, #kg
    m_dot=287, #kg/s
    v_exhaust=3413, #m/s
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

#v2.plot_2d_traj(rel_to_launch_point=False,output_img=True,output_gif=True,gif_frame_count=150)
#mm3.plot_fpa_vs_time(output_img=True)







