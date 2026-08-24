# -*- coding: utf-8 -*-
"""
solve_ode_baseball_traj.py

Description:
    This calculates the 2D trajectory of a baseball using numerical integration.
    
Created by:
    Bradley Canty, 2026/06/24
"""

from numerical_tools import propagate_ode
from matplotlib import pyplot as plt
import numpy as np
import math
import weakref
import os

deg_symbol = u'\N{DEGREE SIGN}'
img_dir = os.getcwd() + r'\Images'
#If the Image directory does not exist, create it
if not os.path.exists(img_dir):
    os.makedirs(img_dir)

class BaseballTraj():
    
    # Class variable to hold references to all active instances
    _instances = weakref.WeakSet()
    
    #BASEBALL PROPERTIES
    Cd = .3 #drag coefficient
    A = math.pi*0.0366**2 #projected frontal area of projectile
    rho_air = 1.225 #kg/m^3, air density
    m = 0.145 #kg, mass of baseball
    g = 9.81 #m/s^2, gravity
    
    def __init__(self,name,v_mag0,launch_angle_deg,final_time_sec,dt,solver):
        BaseballTraj._instances.add(self) #add instance to registry
        self.name = name
        self.v_mag0 = v_mag0
        self.launch_angle_deg = launch_angle_deg
        self.solver = solver

        rx0 = 0 #m, its sensible to keep this at zero
        ry0 = 1 #m, initialize to 1 meter to simulate a person standing with a bat
        vx0 = v_mag0 * np.cos(np.radians(launch_angle_deg))
        vy0 = v_mag0 * np.sin(np.radians(launch_angle_deg))
        self.tspan = [0,final_time_sec]
        self.dt = dt
        state0 = (rx0,ry0,vx0,vy0)
        
        self.times, states = propagate_ode(self.sys_of_ODEs, state0, self.tspan, self.dt, method=self.solver)
        self.rx = states[:,0]
        self.ry = states[:,1]
        self.vx = states[:,2]
        self.vy = states[:,3]
        
        #Check if y value goes below 0, truncate the states at the associated index
        for idx,height in enumerate(self.ry):
            if height < 0:
                self.rx = self.rx[0:idx]
                self.ry = self.ry[0:idx]
                self.vx = self.vx[0:idx]
                self.vy = self.vy[0:idx]
                self.times = self.times[0:idx]
                break

    def sys_of_ODEs(self,t,state):
        '''
        Establish system of ODEs to solve, i.e., the time derivative of the
        state vector
        
        Parameters:
        -----------
        t : float array
            array of times at which to evaluate the function
        state:
            state vector given by (x, y, vx, vy)
        
        Returns:
        -------
        float array
            state vector derivative given by (vx,vy,ax,ay)  
        '''
        x,y,vx,vy = state
        v = math.sqrt(vx ** 2 + vy ** 2) #velocity magnitude
        k_Newton = .5 * BaseballTraj.rho_air * BaseballTraj.Cd * BaseballTraj.A #Newton drag parameter
        zeta = k_Newton / BaseballTraj.m #specific Newton drag parameter
        
        ax = -zeta * v * vx
        ay = -BaseballTraj.g - zeta * v * vy
        
        state_dot = np.array([vx, vy, ax, ay])
        
        return state_dot

    def plot(self):
        plt.plot(self.rx,self.ry,'r',zorder=0)
        plt.title(f'{self.name} Trajectory')
        plt.xlabel('down range [m]')
        plt.ylabel('altitude [m]')
        plt.grid(True)
        plt.axis('equal')
        
        #Plot times at 1 second intervals along the the trajectory
        for i in range(self.tspan[0],len(self.times),int(1/self.dt)):
            plt.scatter(self.rx[i],self.ry[i],color='black',marker='.',s=100,zorder=1,)
            plt.annotate(
                f't={self.times[i]}[s]',
                (self.rx[i], self.ry[i]),
                textcoords="offset points", # How to position the text
                xytext=(0, 10),             # Distance from text to point (x,y) in points
                ha='center'                 # Horizontal alignment: left, right, or center
            )
            
        plt.show()
    
    @classmethod
    def get_all_instances(cls):
        return list(cls._instances)
    
    @classmethod
    def plot_all(cls):
        fig, ax = plt.subplots()
        for traj in list(cls._instances):
            ax.plot(traj.rx,traj.ry,zorder=2,label=f'$v_{{mag,0}}$={round(traj.v_mag0,2)}[m/s]\nlaunch angle={traj.launch_angle_deg}{deg_symbol}')
            
            #Plot times at 1 second intervals along the the trajectory
            for i in range(traj.tspan[0],len(traj.times),int(1/traj.dt)):
                ax.scatter(traj.rx[i],traj.ry[i],color='black',marker='.',s=40,zorder=3)
                ax.annotate(
                    f't={traj.times[i]}[s]',
                    (traj.rx[i], traj.ry[i]),
                    textcoords="offset points",
                    xytext=(0, 5),
                    size=6,
                    ha='center'
                )
                
        ax.set_title('Baseball Trajectories')
        ax.set_xlabel('down range [m]')
        ax.set_ylabel('altitude [m]')
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        ax.grid(True,zorder=0)
        ax.axis('scaled')
        plt.tight_layout()
        plt.savefig(img_dir + '/baseball_trajectories.png', dpi = 600)
        plt.show()

if __name__ == '__main__':
    batted_ball_speed_milesPerHr = 100
    batted_ball_speed_metersPerSec = batted_ball_speed_milesPerHr/2.237
    launch_angles = [20,40,60,80] #deg
    baseballs = []
    for idx,launch_angle in enumerate(launch_angles,1):
        baseball = BaseballTraj(f'Baseball {idx}',batted_ball_speed_metersPerSec,launch_angle,final_time_sec=6,dt=0.1,solver='euler')
        baseballs.append(baseball)
    BaseballTraj.plot_all()




