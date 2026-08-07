# -*- coding: utf-8 -*-
"""
rocket_traj_sim_ex8.py

Description:
    Continuing from exercise 8, extend the simulation to two dimensions. 
    Additionally, implement a pitch program such that the rocket begins to pitch over 
    at an altitude of 30 meters and such that the flight path angle at burnout is
    varied by the user: this initiates the gravity turn.
    
    Plot...
    - 1) the rocket trajectory (i.e. altitude versus downrange distance)
    - 2) vertical distance vs horizontal distance from launch site across limb of earth
    - 3) flight path angle vs time
    - 4) position, velocity, and acceleration vs time
    - 5) dynamic pressure vs time
    - 6) mass vs time
    - 7) Mach number and drag coefficient vs time
    - 8) a GIF of the trajectory

Created by: Bradley Canty, 2026/07/08
"""
import os
import io
import imageio as iio
from PIL import Image
import numpy as np
from symbols import deg_symbol
from matplotlib import pyplot as plt
from numerical_tools import (propagate_traj,
                             get_altitude_km,
                             get_range_km,
                             get_radial_vel,
                             get_tangential_vel,
                             get_radial_accel,
                             get_tangential_accel)
from rocket_constants import m_empty, m_propellant,t_burnout,pitch_programs
from physical_constants import r_earth

img_dir = os.getcwd() + r'\Images'
#If the Image directory does not exist, create it
if not os.path.exists(img_dir):
    os.makedirs(img_dir)

#Initial conditions
lon0_deg = 90 #launch longitude
#lat0_deg = 0 #launch latitude NOT VALID IN 2D SIMULATION
alt0_m = 0 #launch altitude
rx0 = (r_earth+alt0_m)*np.cos(lon0_deg * np.pi/180)
ry0 = (r_earth+alt0_m)*np.sin(lon0_deg * np.pi/180)
r0 = np.array([rx0,ry0,0]) #initial positions
v0 = np.array([0,0,0]) #initial velocities
m0 = m_empty + m_propellant
state0 = [r0[0],r0[1],v0[0],v0[1],m0] #initial state vector

t = np.linspace(0, 500, 1001)

states,metrics = propagate_traj(state0,t,pitch_programs,method='rk4',get_metrics=True)
rx,ry,vx,vy,m = states
ax_grav,ay_grav,ax_thrust,ay_thrust,ax_drag,ay_drag,Cd,M,fpa_deg,eta_deg,dynamic_pressure = metrics
ax = ax_grav + ax_thrust + ax_drag
ay = ay_grav + ay_thrust + ay_drag
altitude_km = get_altitude_km(rx,ry)
range_km = get_range_km(rx,ry)
fpa_rad = fpa_deg * np.pi/180

r = np.array([rx,ry])
v = np.array([vx,vy])
a = np.array([ax,ay])

v_r = get_radial_vel(r,v)
v_t = get_tangential_vel(r,v)
a_r = get_radial_accel(r,a)
a_t = get_tangential_accel(r,a)

'''
GENERATE PLOTS
'''
#1) the rocket trajectory (i.e. altitude versus downrange distance)
plt.figure()
plt.plot(range_km,altitude_km,'r')
plt.title('V2 Rocket Trajectory (ECI frame)')
plt.xlabel('Downrange [km]')
plt.ylabel('Altitude [km]')
plt.grid(True)
plt.axis('equal')
plt.savefig(img_dir + '/1_downrange_vs_alt.png', dpi = 600)
plt.show()

#2) vertical distance vs horizontal distance from launch site across limb of earth
fpa_at_bo_deg = fpa_deg[round(t_burnout/(t[1]-t[0]))]
t_burnout_text = r'$t_{burnout}$'
limb_of_earth = plt.Circle((0, -r_earth/1000), r_earth/1000, color='b')
fig, ax1 = plt.subplots() # note we must use plt.subplots, not plt.subplot
ax1.add_patch(limb_of_earth)

#rotate the coordinates by longitude to make the plot horizontal
lon0_rad = lon0_deg * np.pi/180
R = np.array([[np.cos(-lon0_rad),-np.sin(-lon0_rad)],[np.sin(-lon0_rad),np.cos(-lon0_rad)]])
r_vec = np.array([rx,ry])
r_vec_new = np.matmul(R,r_vec)
rx_new = r_vec_new[0]
if (r_vec_new[1][-1] < 0):
    ry_new = -r_vec_new[1]
else:
    ry_new = r_vec_new[1]

ax1.plot(ry_new/1000,(rx_new-r_earth)/1000,'r')
ax1.set_xlim(min(ry_new)/1000,(max(ry_new) + .25*(max(ry_new)-min(ry_new)))/1000)
ax1.set_ylim((min(rx_new-r_earth) - .25*(max(rx_new-r_earth)-min(rx_new-r_earth)))/1000,(max(rx_new-r_earth) + .25*(max(rx_new-r_earth)-min(rx_new-r_earth)))/1000)
ax1.set_aspect('equal')
ax1.grid()
ax1.set_title(f"V2 Rocket Trajectory Across Limb of Earth (ECI frame)\n\
Pitch program: FPA = {round(pitch_programs['target_fpa_at_bo_deg'],2)}{deg_symbol} (actual {round(fpa_at_bo_deg,2)}{deg_symbol}) at {t_burnout_text} = {round(t_burnout,2)} sec\n\
Max range = {round(max(range_km),2)} [km], Max altitude = {round(max(altitude_km),2)} [km]", fontsize=12)
ax1.set_xlabel('Horizontal dist. from launch site [km]')
ax1.set_ylabel('Vertical dist. from\nlaunch site [km]')

plt.savefig(img_dir + '/2_traj_across_limb_of_earth.png', dpi = 600)
plt.show()

#3) flight path angle vs time
plt.figure()
plt.plot(t[:len(fpa_deg)],fpa_deg,'r',label='flight path angle [deg]')
plt.plot(t[:len(eta_deg)],eta_deg,'g',label='angle from vel. vec. to thrust vec.')
plt.plot(t_burnout, pitch_programs['target_fpa_at_bo_deg'], 'ko',markersize=10,label='target flight path angle at burnout')
plt.axvline(x=t_burnout, color='b', linestyle='--', linewidth=2)
plt.text(x=t_burnout + 0.2, y=-40, s='burnout time', rotation='vertical', color='b',fontsize=15)
plt.suptitle('Flight Path Angle vs Time')
plt.title(f"Pitch program: FPA = {round(fpa_at_bo_deg,2)}{deg_symbol} at {t_burnout_text} = {round(t_burnout,2)} seconds", fontsize=12, color="gray")
plt.xlabel('Time [sec]')
plt.ylabel(f'Angle [{deg_symbol}]')
plt.grid(True)
plt.legend()
#plt.xlim(0,t_burnout+4)
plt.savefig(img_dir + '/3_fpa_vs_time.png', dpi = 600)
plt.show()

#4) position vs time, velocity vs time, and acceleration vs time
fig, (ax1, ax2, ax3) = plt.subplots(1,3,layout="constrained")
fig.suptitle("V2 Rocket Trajectory Metrics (ECI Frame)", fontsize=16, fontweight="bold", y=1.1)
big_ax = fig.add_subplot(111, frameon=False)
big_ax.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
big_ax.grid(False)
t_burnout_text = r'$t_{burnout}$'
big_ax.set_title(f"Pitch program: FPA = {round(fpa_at_bo_deg,2)}{deg_symbol} at {t_burnout_text} = {round(t_burnout,2)} seconds", fontsize=12, color="gray", pad=30, y=1.58)

ax1.plot(t[:len(range_km)]/60, range_km, label='range')
ax1.plot(t[:len(altitude_km)]/60, altitude_km, label='altitude')
ax1.set_title(f'Max range = {round(max(range_km),2)} [km]\nMax altitude = {round(max(altitude_km),2)} [km]',fontsize=8)
ax1.set_xlabel('time [min]')
ax1.set_ylabel('position [km]')
ax1.legend(
    loc="upper center", 
    bbox_to_anchor=(0.5, -0.15), 
    ncol=1
    )
ax1.grid()

ax2.plot(t[:len(v_t)]/60, v_t/1000, label='tangential vel.')
ax2.plot(t[:len(v_r)]/60, v_r/1000, label='radial vel.')
ax2.set_title(f'|Max tangential vel.| = {round(max(abs(v_t/1000)),2)} [km/s]\n|Max radial vel.| = {round(max(abs(v_r/1000)),2)} [km/s]',fontsize=8)
ax2.set_xlabel('time [min]')
ax2.set_ylabel('velocity [km/s]')
ax2.legend(
    loc="upper center", 
    bbox_to_anchor=(0.5, -0.15), 
    ncol=1
    )
ax2.grid()

ax3.plot(t[:len(a_t)]/60, a_t, label='tangential accel.')
ax3.plot(t[:len(a_r)]/60, a_r, label='radial accel.')
ax3.set_title(f'|Max tangential accel.| = {round(max(abs(a_t)),2)} [$m/s^2$]\n|Max radial accel.| = {round(max(abs(a_r)),2)} [$m/s^2$]',fontsize=8)
ax3.set_xlabel('time [min]')
ax3.set_ylabel(r'acceleration [$m/s^2$]')
ax3.legend(
    loc="upper center", 
    bbox_to_anchor=(0.5, -0.15), 
    ncol=1
    )
ax3.grid()
fig.savefig(img_dir + '/4_pos_vel_accel_vs_time.png', dpi = 600)

#6) mass vs time
plt.figure()
plt.plot(t[:len(m)],m,'r')
plt.title('Mass vs Time')
plt.grid(True)
plt.xlabel('time [sec]')
plt.ylabel('mass [kg]')
plt.savefig(img_dir + '/6_mass_vs_time.png', dpi = 600)
plt.show()

#7) Mach number and Drag Coefficient vs Time
plt.figure()
plt.plot(t[:len(M)], M,'r',label='Mach number')
plt.plot(t[:len(Cd)],Cd,'b',label='drag coefficient')
plt.title('Mach Number and Drag Coefficient vs Time')
plt.legend()
plt.grid(True)
plt.savefig(img_dir + '/7_mach_num_and_drag_coeff_vs_time.png', dpi = 600)
plt.show()


'''
plt.figure()
plt.plot(t[:len(ax_thrust)],ax_thrust,'b',label=r'$a_x$ thrust')
plt.plot(t[:len(ay_thrust)],ay_thrust,'r',label=r'$a_y$ thrust')
plt.title(r'$a_{thrust}$ vs Time for V2 Rocket (g = f(h), $\rho_{air} = f(h)\ [kg/m^3]$)')
plt.xlabel('time [s]')
plt.ylabel(r'acceleration [$m/s^2$]')
plt.grid(True)
plt.legend()

plt.figure()
plt.plot(t[:len(ay_drag)],ay_drag,'r')
plt.title(r'$a_{y,drag}$ vs Time for V2 Rocket (g = f(h), $\rho_{air} = f(h)\ [kg/m^3]$)')
plt.xlabel('time [s]')
plt.ylabel(r'acceleration [$m/s^2$]')
plt.grid(True)

plt.figure()
plt.plot(t[:len(ay_grav)],ay_grav,'r')
plt.title(r'$a_{y,grav}$ vs Time for V2 Rocket (g = f(h), $\rho_{air} = f(h)\ [kg/m^3]$)')
plt.xlabel('time [s]')
plt.ylabel(r'acceleration [$m/s^2$]')
plt.grid(True)

plt.figure()
plt.plot(t[:len(ax)],ax,'r',label='$a_x$')
plt.plot(t[:len(ay)],ay,'b',label='$a_y$')
plt.title(r'acceleration vs Time for V2 Rocket (g = f(h), $\rho_{air} = f(h)\ [kg/m^3]$)')
plt.xlabel('time [s]')
plt.ylabel(r'acceleration [$m/s^2$]')
plt.legend()
plt.grid(True)

plt.figure()
plt.plot(t[:len(vx)],vx,'r',label='$v_x$')
plt.plot(t[:len(vy)],vy,'b',label='$v_y$')
plt.title(r'velocity vs Time (g = f(h), $\rho_{air} = f(h)\ [kg/m^3]$)')
plt.xlabel('time [s]')
plt.ylabel('velocity [m/s]')
plt.legend()
plt.grid(True)

plt.figure()
plt.plot(t[:len(vx)]/60,rx/1000,'r',label='$r_x$')
plt.plot(t[:len(vy)]/60,(ry-r_earth)/1000,'b',label='$r_y$')
plt.title(r'Position vs Time (g = f(h), $\rho_{air} = f(h)\ [kg/m^3]$)')
plt.xlabel('time [min]')
plt.ylabel('position [km]')
plt.legend()
plt.grid(True)

plt.figure()
plt.plot(t[:len(altitude_km)],altitude_km,'r')
plt.title(r'Altitude vs Time (g = f(h), $\rho_{air} = f(h)\ [kg/m^3]$)')
plt.xlabel('time [s]')
plt.ylabel('Altitude [km]]')
plt.grid(True)
'''

#8) Generate a GIF of vertical distance vs horizontal distance from launch site across limb of earth
def plot_traj(i,img_name):
    fpa_at_bo_deg = fpa_deg[round(t_burnout/(t[1]-t[0]))]
    t_burnout_text = r'$t_{burnout}$'
    limb_of_earth = plt.Circle((0, -r_earth/1000), r_earth/1000, color='b')
    fig, ax1 = plt.subplots()
    ax1.add_patch(limb_of_earth)
    
    #rotate the coordinates by longitude to make the plot horizontal
    lon0_rad = lon0_deg * np.pi/180
    R = np.array([[np.cos(-lon0_rad),-np.sin(-lon0_rad)],[np.sin(-lon0_rad),np.cos(-lon0_rad)]])
    r_vec = np.array([rx,ry])
    r_vec_new = np.matmul(R,r_vec)
    rx_new = r_vec_new[0]
    if (r_vec_new[1][-1] < 0):
        ry_new = -r_vec_new[1]
    else:
        ry_new = r_vec_new[1]

    ax1.plot(ry_new[0:i]/1000,(rx_new[0:i]-r_earth)/1000,'r')
    ax1.plot(ry_new[i]/1000,(rx_new[i]-r_earth)/1000,'y.',markersize='10')
    ax1.set_xlim(min(ry_new)/1000,(max(ry_new) + .25*(max(ry_new)-min(ry_new)))/1000)
    ax1.set_ylim((min(rx_new-r_earth) - .25*(max(rx_new-r_earth)-min(rx_new-r_earth)))/1000,(max(rx_new-r_earth) + .25*(max(rx_new-r_earth)-min(rx_new-r_earth)))/1000)
    ax1.set_aspect('equal')
    ax1.grid()
    ax1.set_title(f"V2 Rocket Trajectory Across Limb of Earth (ECI frame)\nPitch program: FPA = {round(pitch_programs['target_fpa_at_bo_deg'],2)}{deg_symbol} (actual {round(fpa_at_bo_deg,2)}{deg_symbol}) at {t_burnout_text} = {round(t_burnout,2)} sec", fontsize=12)
    ax1.set_xlabel('Horizontal dist. from launch site [km]')
    ax1.set_ylabel('Vertical dist. from\nlaunch site [km]')
    
    # Save the plot to an in-memory buffer
    buf = io.BytesIO()
    fig.savefig(buf,format='png',bbox_inches='tight')
    buf.seek(0)
    
    plt.close()
    
    # Append the PIL image to our list of frames
    frames.append(Image.open(buf))

n = 30 #number of images to compose the GIF
count = 0
frames = []
for i in range(0,len(rx),round(len(rx)/(n-2))):
    img_name = f'v2_traj_{count}.png'
    plot_traj(i,img_name)
    count += 1

#generate final image
i = len(rx)-1
img_name = f'v2_traj_{count}.png'
plot_traj(i,img_name)

#Save as GIF using imageio.mimsave (analogous to imwrite)
iio.mimsave('Images/8_v2_traj_gif.gif', frames, format='GIF', 
                duration=100, loop=1) # loop=0 sets an infinite loop

print(f'Trajectory plots found at:\n{img_dir}')
