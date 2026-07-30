# -*- coding: utf-8 -*-
"""
Created on Thu Jul  9 21:22:07 2026

@author: bradl


"""
import io
import imageio as iio
from matplotlib import pyplot as plt
from PIL import Image
from physical_constants import r_earth, OMEGA_EARTH
import numpy as np
import numerical_tools as nt
import sys
from symbols import deg_symbol, eta_symbol, gamma_symbol
import os

img_dir = os.getcwd() + r'\Images'
#If the Image directory does not exist, create it
if not os.path.exists(img_dir):
    os.makedirs(img_dir)

class RocketStage():
    def __init__(self,name,m_empty,m_propellant,m_dot,v_exhaust,stage_diameter_m,coast_duration_sec,t_steer_start_rel,steer_duration_sec,target_fpa_at_bo_deg):
        
        self.name = name
        self.m_empty = m_empty
        self.m_propellant = m_propellant
        self.m_dot = m_dot
        self.v_exhaust = v_exhaust
        self.stage_diameter_m = stage_diameter_m
        self.t_steer_start_rel = t_steer_start_rel #steer start time relative to engine start, with engine starting AFTER the coast duration
        
        self.target_fpa_at_bo_deg = target_fpa_at_bo_deg
        self.coast_duration_sec = coast_duration_sec
        
        #Derived variables
        self.A = np.pi * (self.stage_diameter_m)**2 / 4
        self.engine_burn_duration_sec = self.m_propellant/self.m_dot #sec, burnout time relative to stage's engine start
        
        if steer_duration_sec > self.engine_burn_duration_sec:
            self.steer_duration_sec = self.engine_burn_duration_sec
        elif steer_duration_sec == -1:
            self.steer_duration_sec = self.engine_burn_duration_sec
        else:
            self.steer_duration_sec = steer_duration_sec
            
        #Variables to be computed upon program execution
        self.required_eta_deg = None #temporary, should be None
        self.fpa_at_bo_deg = None #temporary, should be None
        self.t_start_abs = None #absolute start time, i.e. start time relative to t-0
        self.t_bo_abs = None #absolute burnout time
        self.m_other = None #kg, the mass of all other upper stages and payload
        
def null_config():
	return {
        'name'               : 'null',
		'date0'              : '2020-01-01', #initial date
		'et0'                : None,         #initial ephemeris time, i.e., seconds since J2000
		'frame'              : 'J2000',      #standard for earth-centered propagation
		'ode_solver'         : 'rk4',
        'lla0'               : np.array([0,0,0]), #Latitude [degrees], Longitude [degrees], altitude [meters]
        'v0_rel'             : np.array([0,0,0]), #vx [m/s], vy [m/s], vz [m/s]
		'stages'             : [],
        'tspan'              : 0,
        'dt'                 : 0,
        'm_payload'          : 0,
        'launch_azimuth_deg' : 0,
		'output_dir'         : f'{img_dir}',
	}

class Rocket():
    def __init__(self,config):
        self.config = null_config()
        for key in config.keys():
            self.config[key] = config[key]
        
        if self.config['lla0'][0] !=0:
            print(f"ERROR: for this 2D rocket simulation, the launch site latitude should be set to 0{deg_symbol}! Right now, the {self.config['name']} launch site latitude is set to {self.config['lla0'][0]}{deg_symbol}. Please set {self.config['name']} launch site latitude to 0{deg_symbol}!")
            sys.exit(1) #end program with exit code 1
        
        if abs(self.config['launch_azimuth_deg']) != 90:
            print(f"ERROR: for this 2D rocket simulation, the launch azimuth should be either 90 or -90 degrees! Right now, launch_azimuth_deg = {self.config['launch_azimuth_deg']}")
            sys.exit(1) #end program with exit code 1
            
        #If the coast time of the first stage is set to a non-zero value, warn the user to check the inputs and end the program
        if self.config['stages'][0].coast_duration_sec != 0:
            print(f"ERROR: the {self.config['name']} first stage coast time is {self.config['stages'][0].coast_duration_sec} sec while it should be 0! Please review your stage configurations.")
            sys.exit(1) #end program with exit code 1
        
        self.lat0_deg = self.config['lla0'][0]
        self.lon0_deg = self.config['lla0'][1]
        self.alt0_m = self.config['lla0'][2]
        
        self.rx0 = (r_earth+self.alt0_m)*np.cos(self.lon0_deg * np.pi/180)
        self.ry0 = (r_earth+self.alt0_m)*np.sin(self.lon0_deg * np.pi/180)
        
        self.r0 = np.array([self.rx0,self.ry0,0]) #initial positions
        self.v0_rel = np.array([0,0,0]) #initial velocities relative to launch pad
        omega_earth = np.array([0,0,OMEGA_EARTH])
        self.v0 = self.v0_rel + np.cross(omega_earth,self.r0) #initial velocity vector, accounting for rotating atmosphere
        
        self.solve_stage_times()
        self.solve_stage_mass_other()
        
        self.m0 = self.config['stages'][0].m_empty + self.config['stages'][0].m_propellant + self.config['stages'][0].m_other
        self.state0 = [self.r0[0],self.r0[1],self.v0[0],self.v0[1],self.m0] #initial state vector
        
        states,metrics = nt.propagate_traj(self.state0,
                                           self.config['tspan'],
                                           self.config['dt'],
                                           self.config['stages'],
                                           self.config['launch_azimuth_deg'],
                                           self.config['payload_separates'],
                                           self.config['ode_solver'],
                                           get_metrics=True)
        
        self.rx = states[0]
        self.ry = states[1]
        self.vx = states[2]
        self.vy = states[3]
        self.m =  states[4]

        self.ax_grav = metrics[0]
        self.ay_grav = metrics[1]
        self.ax_thrust = metrics[2]
        self.ay_thrust = metrics[3]
        self.ax_drag = metrics[4]
        self.ay_drag = metrics[5]
        self.Cd = metrics[6]
        self.M = metrics[7]
        self.fpa_deg = metrics[8]
        self.eta_deg = metrics[9]
        self.dynamic_pressure = metrics[10]
        
        self.ax = self.ax_grav + self.ax_thrust + self.ax_drag
        self.ay = self.ay_grav + self.ay_thrust + self.ay_drag
        self.fpa_rad = self.fpa_deg * np.pi/180
        
        self.r = np.array([self.rx,self.ry])
        self.v = np.array([self.vx,self.vy])
        self.a = np.array([self.ax,self.ay])
        
        self.altitude_km = nt.get_altitude_km(self.rx,self.ry)
        self.range_km = nt.get_range_km(self.rx,self.ry)
        self.v_rad = nt.get_radial_vel(self.r,self.v)
        self.v_tan = nt.get_tangential_vel(self.r,self.v)
        self.a_rad = nt.get_radial_accel(self.r,self.a)
        self.a_tan = nt.get_tangential_accel(self.r,self.a)
        
        self.stage_count = len(self.config['stages'])
        self.t_bo_final = round(self.config['stages'][self.stage_count-1].t_bo_abs)
        idx_t_bo_final = round(self.t_bo_final/self.config['dt'])
        self.v_tan_bo = self.v_tan[idx_t_bo_final] #tangential velocity at final stage burnout
        
        self.times = np.arange(0,self.config['tspan']+self.config['dt']/2,self.config['dt'])[:len(states[0])]
        self.output_msg_displayed_once = False
        
    def solve_stage_times(self):
        '''
        This solve for the stage start time using the sum of the coast and engine burn durations of previous stages
        '''
        for i,stage in enumerate(self.config['stages']):
            #print(f'i = {i}, stage.name = {stage.name}')
            stage.t_start_abs = 0 
            stage.t_bo_abs = 0
            #Sum the previous coast and engine burn durations
            if i > 0:
                stage.t_start_abs += self.config['stages'][i-1].t_start_abs + self.config['stages'][i-1].coast_duration_sec + self.config['stages'][i-1].engine_burn_duration_sec
                
            stage.t_bo_abs += stage.t_start_abs + self.config['stages'][i].engine_burn_duration_sec
            
    def solve_stage_mass_other(self):
        '''
        This solves for the "other mass" which is the sum of the mass of all the stages above the current stage 
        and the payload mass. Thus, in the simulation the total mass to be used upon staging is:
        m_empty + m_propellant + m_other
        '''
        m_payload = self.config['m_payload']
        stages = self.config['stages']
        for i in range(len(stages)):
            stages[i].m_other = m_payload
            for j in range(i+1,len(self.config['stages'])):
                stages[i].m_other += stages[j].m_empty + stages[j].m_propellant
            #print(f'Stage {i+1} m_other = {stages[i].m_other}')
                
        
    def plot_alt_vs_range(self,output_img=False):
        plt.figure()
        plt.plot(self.range_km,self.altitude_km,'r')
        plt.ylim(0, None) 
        plt.title(f'{self.config['name']} Rocket Trajectory (ECI frame)')
        plt.xlabel('Range [km]')
        plt.ylabel('Altitude [km]')
        plt.grid(True)
        plt.axis('scaled')
        
        if output_img == True:
            plt.savefig(img_dir + f'/{self.config['name']}_alt_vs_range.png', dpi = 600)
        plt.show()
        
        if output_img and self.output_msg_displayed_once == False:
            print(f'{self.config['name']} media found at found at:\n{self.config['output_dir']}\n')
            self.output_msg_displayed_once = True

    def plot_fpa_vs_time(self,output_img=False):
        plt.figure()
        plt.plot(self.times,self.fpa_deg,'r',label=f'flight path angle, {gamma_symbol} [{deg_symbol}]')
        plt.plot(self.times,self.eta_deg,'g',label=f'steering angle, {eta_symbol} [{deg_symbol}]')
    
        plt.title(f'{self.config['name']} Pitch Program')
        plt.xlabel('Time [s]')
        plt.ylabel(f'Angle [{deg_symbol}]')
        plt.grid(True)
    
        color_dict = {0:'k',1:'c',2:'m'}
        for i,stage in enumerate(self.config['stages']):
            stage_bo_time = stage.t_start_abs + stage.coast_duration_sec + stage.engine_burn_duration_sec
            plt.scatter(stage_bo_time, stage.target_fpa_at_bo_deg, color = color_dict[i],marker='o',s=50,label=f'Stage {i+1} target FPA at $t_{{bo}}$={round(stage_bo_time,1)} [s]',zorder=2)
            x_line_pos = stage.t_start_abs + stage.coast_duration_sec + stage.engine_burn_duration_sec
            plt.axvline(x=x_line_pos, color=color_dict[i], linestyle='--', linewidth=1)
        plt.legend()
        
        if output_img == True:
            plt.savefig(img_dir + f'/{self.config['name']}_fpa_vs_time.png', dpi = 600)
        plt.show()
        
        if output_img and self.output_msg_displayed_once == False:
            print(f'{self.config['name']} media found at found at:\n{self.config['output_dir']}\n')
            self.output_msg_displayed_once = True
    
    def plot_alt_vs_time(self,output_img=False):
        #Altitude vs time
        plt.figure()
        plt.plot(self.times,self.altitude_km,'r')
        plt.title(f'{self.config['name']} Altitude vs Time')
        plt.xlabel('Time [sec]')
        plt.ylabel('Altitude [km]')
        plt.grid(True)
        
        if output_img == True:
            plt.savefig(img_dir + f'/{self.config['name']}_alt_vs_time.png', dpi = 600)
        plt.show()
        
        if output_img and self.output_msg_displayed_once == False:
            print(f'{self.config['name']} media found at found at:\n{self.config['output_dir']}\n')
            self.output_msg_displayed_once = True
    
    def plot_pos_vel_accel_vs_time(self,output_img=False):
        #4) position vs time, velocity vs time, and acceleration vs time
        fig, (ax1, ax2, ax3) = plt.subplots(1,3)
        fig.suptitle(f'{self.config['name']} Trajectory Metrics (ECI Frame)', fontsize=16, fontweight="bold",y=1.02)
        
        big_ax = fig.add_subplot(111, frameon=False)
        big_ax.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
        big_ax.grid(False)
    
        ax1.plot(self.times/60, self.range_km, label='range')
        ax1.plot(self.times/60, self.altitude_km, label='altitude')
        ax1.set_title(f'Max range = {round(max(self.range_km),1)} [km]\nMax altitude = {round(max(self.altitude_km),1)} [km]',fontsize=7)
        ax1.set_xlabel('time [min]')
        ax1.set_ylabel('position [km]')
        ax1.legend(
            loc="upper center", 
            bbox_to_anchor=(0.5, -0.15), 
            ncol=1
            )
        ax1.grid()
    
        ax2.plot(self.times/60, self.v_tan/1000, label='tangential vel.')
        ax2.plot(self.times/60, self.v_rad/1000, label='radial vel.')
        ax2.set_title(f'|Max tangential vel| = {round(max(abs(self.v_tan/1000)),1)} [km/s]\n|Max radial vel| = {round(max(abs(self.v_rad/1000)),1)} [km/s]',fontsize=7)
        ax2.set_xlabel('time [min]')
        ax2.set_ylabel('velocity [km/s]')
        ax2.legend(
            loc="upper center", 
            bbox_to_anchor=(0.5, -0.15), 
            ncol=1
            )
        ax2.grid()
    
        ax3.plot(self.times/60, self.a_tan, label='tangential accel.')
        ax3.plot(self.times/60, self.a_rad, label='radial accel.')
        ax3.set_title(f'|Max tangential accel| = {round(max(abs(self.a_tan)),1)} [$m/s^2$]\n|Max radial accel| = {round(max(abs(self.a_rad)),1)} [$m/s^2$]',fontsize=7)
        ax3.set_xlabel('time [min]')
        ax3.set_ylabel(r'acceleration [$m/s^2$]')
        ax3.legend(
            loc="upper center", 
            bbox_to_anchor=(0.5, -0.15), 
            ncol=1
            )
        ax3.grid()
        
        plt.subplots_adjust(hspace=0.5, wspace=0.6)
        
        if output_img == True:
            plt.savefig(img_dir + f'/{self.config['name']}_pos_vel_accel_vs_time.png', bbox_inches='tight', dpi = 600)
        plt.show()
        
        if output_img and self.output_msg_displayed_once == False:
            print(f'{self.config['name']} media found at found at:\n{self.config['output_dir']}\n')
            self.output_msg_displayed_once = True
    
    def plot_dyn_pres_vs_time(self,output_img=False):
        #5) dynamic pressure vs time
        plt.figure()
        plt.plot(self.times,self.dynamic_pressure/1000,'r')
        plt.title(f'{self.config['name']} Dynamic Pressure vs Time')
        plt.grid(True)
        plt.xlabel('time [s]')
        plt.ylabel(r'dynamic pressure [KPa]')
        
        if output_img == True:
            plt.savefig(img_dir + f'/{self.config['name']}_dyn_pres_vs_time.png', dpi = 600)
        plt.show()
        
        if output_img and self.output_msg_displayed_once == False:
            print(f'{self.config['name']} media found at found at:\n{self.config['output_dir']}\n')
            self.output_msg_displayed_once = True
    
    def plot_mass_vs_time(self,output_img=False):
        plt.figure()
        plt.plot(self.times,self.m,'r')
        plt.title(f'{self.config['name']} Mass vs Time')
        plt.grid(True)
        plt.xlabel('time [s]')
        plt.ylabel('mass [kg]')
        
        if output_img == True:
            plt.savefig(img_dir + f'/{self.config['name']}_mass_vs_time.png', dpi = 600)
        plt.show()
        
        if output_img and self.output_msg_displayed_once == False:
            print(f'{self.config['name']} media found at found at:\n{self.config['output_dir']}\n')
            self.output_msg_displayed_once = True
    
    def plot_mach_cd_vs_time(self,output_img=False):
        plt.figure()
        plt.plot(self.times, self.M,'r',label='Mach number')
        plt.plot(self.times,self.Cd,'b',label='drag coefficient')
        plt.xlabel('time [s]')
        plt.title(f'{self.config['name']} Mach Number and Drag Coefficient vs Time')
        plt.legend()
        plt.grid(True)
        
        if output_img == True:
            plt.savefig(img_dir + f'/{self.config['name']}_mach_cd_vs_time.png', dpi = 600)
        plt.show()
        
        if output_img and self.output_msg_displayed_once == False:
            print(f'{self.config['name']} media found at found at:\n{self.config['output_dir']}\n')
            self.output_msg_displayed_once = True

    def plot_2d_traj(self,rel_to_launch_point=True,output_img=False,output_gif=False,gif_frame_count=30,ms_per_frame=100):
        #Generate plot of vertical distance vs horizontal distance from ECI launch point across limb of earth
        
        if rel_to_launch_point == True:
            limb_of_earth = plt.Circle((0, -r_earth/1000), r_earth/1000, color='b',zorder=0)
            fig, ax1 = plt.subplots()
            ax1.add_patch(limb_of_earth)
    
            #rotate the coordinates by longitude to make the plot horizontal
            lon0_rad = self.lon0_deg * np.pi/180
            R = np.array([[np.cos(-lon0_rad),-np.sin(-lon0_rad)],[np.sin(-lon0_rad),np.cos(-lon0_rad)]])
            r_vec = np.array([self.rx,self.ry])
            r_vec_new = np.matmul(R,r_vec)
            rx_new = r_vec_new[0]
            if (r_vec_new[1][-1] < 0):
                ry_new = -r_vec_new[1]
            else:
                ry_new = r_vec_new[1]
    
            ax1.plot(ry_new/1000,(rx_new-r_earth)/1000,'r', zorder=2)
            ax1.scatter(ry_new[0]/1000,(rx_new[0]-r_earth)/1000,c='cyan',marker='.',s=100,label='launch point', zorder=3)
            if self.altitude_km[-1] < 2:
                ax1.scatter(ry_new[-1]/1000,(rx_new[-1]-r_earth)/1000,c='yellow',marker='.',s=100,label='impact point', zorder=3)
            elif self.altitude_km[-1] > 2:
                ax1.scatter(ry_new[-1]/1000,(rx_new[-1]-r_earth)/1000,c='orange',marker='.',s=100,label='final point', zorder=3)
            ax1.set_xlim((min(ry_new) - .25*(max(ry_new)-min(ry_new)))/1000,(max(ry_new) + .25*(max(ry_new)-min(ry_new)))/1000)
            ax1.set_ylim((min(rx_new-r_earth) - .25*(max(rx_new-r_earth)-min(rx_new-r_earth)))/1000,(max(rx_new-r_earth) + .25*(max(rx_new-r_earth)-min(rx_new-r_earth)))/1000)
            ax1.set_xlabel('Horizontal dist. from ECI launch point [km]')
            ax1.set_ylabel('Vertical dist. from\nECI launch point [km]')
        else:
            limb_of_earth = plt.Circle((0, 0), r_earth/1000, color='b')
            fig, ax1 = plt.subplots()
            ax1.add_patch(limb_of_earth)
            ax1.plot(self.rx/1000,self.ry/1000,'r', zorder=2)
            ax1.scatter(self.rx[0]/1000,self.ry[0]/1000,c='cyan',marker='.',s=100,label='launch point', zorder=3)
            if self.altitude_km[-1] < 2:
                ax1.scatter(self.rx[-1]/1000,self.ry[-1]/1000,c='yellow',marker='.',s=100,label='impact point', zorder=3)
            if self.altitude_km[-1] > 2:
                ax1.scatter(self.rx[-1]/1000,self.ry[-1]/1000,c='orange',marker='.',s=100,label='final point', zorder=3)
            ax1.set_xlim((min(self.rx) - .25*(max(self.rx)-min(self.rx)))/1000,(max(self.rx) + .25*(max(self.rx)-min(self.rx)))/1000)
            ax1.set_ylim((min(self.ry) - .25*(max(self.ry)-min(self.ry)))/1000,(max(self.ry) + .25*(max(self.ry)-min(self.ry)))/1000)
            ax1.set_xlabel('ECI x-coordinate [km]')
            ax1.set_ylabel('ECI y-coordinate [km]')
        
        ax1.legend(fontsize=8)
        ax1.set_aspect('equal')
        ax1.grid(zorder=1)    
        ax1.set_title(f"{self.config['name']} Trajectory Across Limb of Earth (ECI frame)\nMax range = {round(max(self.range_km),1)} [km], Max alt = {round(max(self.altitude_km),1)} [km]", fontsize=12)
        
        if output_img == True:
            plt.savefig(img_dir + f'/{self.config['name']}_2d_traj.png', dpi = 600)
        plt.show()

        def get_2d_traj_frame(i):
            if rel_to_launch_point == True:
                limb_of_earth = plt.Circle((0, -r_earth/1000), r_earth/1000, color='b',zorder=0)
                fig, ax1 = plt.subplots()
                ax1.add_patch(limb_of_earth)
                
                #rotate the coordinates by longitude to make the plot horizontal
                lon0_rad = self.lon0_deg * np.pi/180
                R = np.array([[np.cos(-lon0_rad),-np.sin(-lon0_rad)],[np.sin(-lon0_rad),np.cos(-lon0_rad)]])
                r_vec = np.array([self.rx,self.ry])
                r_vec_new = np.matmul(R,r_vec)
                rx_new = r_vec_new[0]
                if (r_vec_new[1][-1] < 0):
                    ry_new = -r_vec_new[1]
                else:
                    ry_new = r_vec_new[1]
        
                ax1.plot(ry_new[0:i]/1000,(rx_new[0:i]-r_earth)/1000,'r', zorder=2)
                ax1.scatter(ry_new[0]/1000,(rx_new[0]-r_earth)/1000,c='cyan',marker='.',s=100,label='launch point', zorder=3)
                if i == len(self.rx)-1 and self.altitude_km[i] < 2:
                    ax1.scatter(ry_new[-1]/1000,(rx_new[-1]-r_earth)/1000,color='yellow',marker='.',s=100,label='impact point', zorder=3)
                elif i == len(self.rx)-1 and self.altitude_km[i] > 2:
                    ax1.scatter(ry_new[i]/1000,(rx_new[i]-r_earth)/1000,color='orange',marker='.',s=100,label='final point', zorder=3)
                else:
                    ax1.scatter(ry_new[i]/1000,(rx_new[i]-r_earth)/1000,color='orange',marker='.',s=100, zorder=3)
                    
                ax1.set_xlim((min(ry_new) - .25*(max(ry_new)-min(ry_new)))/1000,(max(ry_new) + .25*(max(ry_new)-min(ry_new)))/1000)
                ax1.set_ylim((min(rx_new-r_earth) - .25*(max(rx_new-r_earth)-min(rx_new-r_earth)))/1000,(max(rx_new-r_earth) + .25*(max(rx_new-r_earth)-min(rx_new-r_earth)))/1000)
                ax1.set_xlabel('Horizontal dist. from ECI launch point [km]')
                ax1.set_ylabel('Vertical dist. from\n ECI launch point [km]')
            else:
                limb_of_earth = plt.Circle((0, 0), r_earth/1000, color='b',zorder=0)
                fig, ax1 = plt.subplots()
                ax1.add_patch(limb_of_earth)
                
                ax1.plot(self.rx[0:i]/1000,self.ry[0:i]/1000,'r', zorder=2)
                ax1.scatter(self.rx[0]/1000,self.ry[0]/1000,c='cyan',marker='.',s=100,label='launch point', zorder=3)
                if i == len(self.rx)-1 and self.altitude_km[i] < 2:
                    ax1.scatter(self.rx[-1]/1000,self.ry[-1]/1000,c='yellow',marker='.',s=100,label='impact point', zorder=3)
                elif i == len(self.rx)-1 and self.altitude_km[i] > 2:
                    ax1.scatter(self.rx[i]/1000,self.ry[i]/1000,color='orange',marker='.',s=100,label='final point', zorder=3)
                else:
                    ax1.scatter(self.rx[i]/1000,self.ry[i]/1000,color='orange',marker='.',s=100, zorder=3)
                  
                ax1.set_xlim((min(self.rx) - .25*(max(self.rx)-min(self.rx)))/1000,(max(self.rx) + .25*(max(self.rx)-min(self.rx)))/1000)
                ax1.set_ylim((min(self.ry) - .25*(max(self.ry)-min(self.ry)))/1000,(max(self.ry) + .25*(max(self.ry)-min(self.ry)))/1000)
                ax1.set_xlabel('ECI x-coordinate [km]')
                ax1.set_ylabel('ECI y-coordinate [km]')
            
            ax1.legend(fontsize=8)
            ax1.set_aspect('equal')
            ax1.grid(zorder=1)
            whitespace = '  '
            ax1.set_title(f"{whitespace}{self.config['name']} Trajectory Across Limb of Earth (ECI frame){whitespace}\nAt t+{self.times[i]}[s]: alt={round(self.altitude_km[i],1)}[km], vel={round(np.sqrt(self.v[0,i]**2 + self.v[1,i]**2),1)}[m/s]", fontsize=12)
            
            # Save the plot to an in-memory buffer
            buf = io.BytesIO()
            fig.savefig(buf,format='png',bbox_inches='tight')
            buf.seek(0)
            plt.close()
            
            if Image.open(buf).size[0] == 599:
                print(f'The frame width is {Image.open(buf).size[0]}')
                pass
            
            return buf
        
        if output_gif == True:
            count = 0
            frames = []
            for i in range(0,len(self.rx),round(len(self.rx)/(gif_frame_count-2))):
                buf = get_2d_traj_frame(i)
                frames.append(Image.open(buf))
                count += 1
        
            #generate final frame
            i = len(self.rx)-1
            buf = get_2d_traj_frame(i)
            frames.append(Image.open(buf))
            
            #Save as GIF using imageio.mimsave (analogous to imwrite)
            gif_file_name = f'\{self.config['name']}_2d_traj_gif.gif'
            gif_file_path = img_dir + gif_file_name
            iio.mimsave(gif_file_path, frames, format='GIF', 
                            duration=ms_per_frame, loop=1) # loop=0 sets an infinite loop
            
        if (output_img or output_gif) and self.output_msg_displayed_once == False:
            print(f'{self.config['name']} media found at found at:\n{self.config['output_dir']}\n')
            self.output_msg_displayed_once = True

    def generate_all_plots(self,rel_to_launch_point=False,output_img=False,output_gif=False,gif_frame_count=50):
        self.plot_alt_vs_time(output_img=True)
        self.plot_fpa_vs_time(output_img=True)
        self.plot_pos_vel_accel_vs_time(output_img=True)
        self.plot_dyn_pres_vs_time(output_img=True)
        self.plot_mass_vs_time(output_img=True)
        self.plot_mach_cd_vs_time(output_img=True)
        self.plot_2d_traj(rel_to_launch_point,output_img,output_gif,gif_frame_count)

