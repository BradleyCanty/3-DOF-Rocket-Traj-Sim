# 3-DOF Rocket Trajectory Simulation Tutorial [WIP]
## Summary
This is a step-by-step tutorial for creating a 3-degree-of-freedom rocket trajectory simulation in Python. It begins by comparing a numerical simulation of 1D rocket motion with the Tsiolkovsky rocket equation, and concludes with a desktop application that takes rocket properties as input (launch location, launch azimuth, launch time, mass parameters per stage, engine parameters per stage, desired burnout flight path angle per stage, payload mass, etc.) and outputs... 
1) plots showing the rocket's trajectory metrics
2) an interactive 3D visualization of the rocket's trajectory over an oblate, rotating Earth

## Learning Objectives
Upon completion you should be able to...
1) formulate the equations of motion for a rocket
2) solve the equations of motion to get position and velocity at every time step along its trajectory
3) understand the reference frames and coordinate systems used in astrodynamics
4) compute the six classical orbital elements from the rocket's state vector at burnout
5) understand how Earth's oblateness affects the rocket's motion

## Simulation Architecture
[COMPLETE THIS SECTION]

## Usage and Installation
All exercises are implemented in python and are available in the GitHub repo.
Each exercise has an associated main file (rocket_traj_sim_ex.py):
* the header 'Description' section contains instructions for completing the exercise
* the body contains the solution to the exercise; as a student it is incumbent on you to attempt a solution on your own: if you get stuck, simply refer to the working code

Install the dependencies using the 'requirements.txt' file:
1) Open your terminal or command prompt
2) Navigate to the directory containing the file:\
   cd path/to/project
5) Activate a new virtual environment for this project (strongly recommended to avoid messing up your global Python installation)
6) Run the installation command:\
   pip install -r requirements.txt

## Outline of Exercises
1) Simulate the 1-D motion of V2 rocket in absence of outside forces using Euler numerical integration method. Compare results with those of the Tsiolkovsky rocket equation.
2) Modify exercise 1 to apply constant gravity and a thrust cutoff when propellant is expended
3) Apply altitude-dependent gravity, and solve the system of ODEs using midpoint method and 4th order Runge Kutta method
4) Add in the effects of drag using a constant drag coefficient Cd = 0.125
5) Compute air density as a function of altitude using the isothermal barometric formula. Then, simulate the trajectory of the rocket using this altitude-varying air density formula.
6) Fit speed of sound versus altitude data with piecewise continuous lines. Additionally, fit drag coefficient vs Mach number data (derived from V2 wind tunnel tests) with a natural cubic spline.
7) Implement the Mach number-varying drag coefficient in the 1D rocket trajectory simulation using the curve fits from the previous exercise.
8) Extend the rocket trajectory simulation from 1D to 2D, and apply an initial offset from vertical to initiate a gravity turn. Additionally, make plots of...
   * altitude versus downrange distance
   * flight path angle vs time
   * velocity vs time
   * acceleration vs time
   * dynamic pressure vs time
   * mass vs time
   * trajectory and limb of earth
 9) Starting from the previous exercise, replace the initial offset from vertical with a vertical lift-off followed by a pitch over (i.e. implement a "pitch program"). Additionally, make the same plots as before.
10) Assuming an Earth Centered Inertial frame in the simulation we've built so far, apply the velocity of the rocket relative to the Earth's rotating atmosphere. This involves reformulating the flight path angle calculation.
11) Organize the code to instantiate a 'rocket' object for each simulation. Organize the plotting functionality such that plots are generated using method calls on the rocket objects. Additionally, extend the code to handle rockets with multiple stages. Finally, make a plot showing the trajectory of a multi-stage rocket.
12) Split the rocket trajectory flight into phases and show the phases in the plots. This will make it easier to discern what's going on:\
    Stage 1:
    1) vertical flight
    2) pitch over
    3) gravity turn
   
    For j in remaining stages:
    1) jth stage burn
    2) jth engine cutoff
    
    Update all plots to show the trajectory colored by the flight phases, and label the phases using a legend.

13) Refactor the atmospheric model to use the U.S. Standard Atmosphere 1976 specification, then calculate pressure as a function of height, p(h). Update the density, temperature, and speed of sound calculations to use this more accurate atmospheric model. Additionally, update the rocket thrust equation to account for the change in atmospheric pressure, that is, account for the "pressure thrust".
The current rocket thrust equation only uses mass flow rate, m_dot, and exhaust velocity, v_exhaust:\
T = m_dot * v_exhaust\
A more physically-accurate rocket thrust equation introduces a term to account for the change in atmospheric pressure with altitude:\
T = m_dot * v_exhaust + (p_exhaust - p(h)) * A_exit\
where\
p_exhaust = engine exhaust gas pressure at the engine exit plane (i.e. at the end of the nozzle)
p(h)      = air pressure as a function of altitude h
A_exit    = area at the engine exit plane (i.e. at the end of the nozzle)\
Update the RocketStage class to have the p_exhaust and A_exit properties required for computing thrust, then update the thrust terms in the get_state_dot() function in the numerical_tools.py file

14) Implement the following coordinate transformations:
    * LLA to ECI
    * ECI to LLA
    * LLA to ECEF
    * ECEF to LLA
    
    Implement the launch location by starting in the Latitude, Longitude, Altitude (LLA) frame, converting it to the Earth Centered, Earth Fixed (ECEF) frame, and then converting it to the Earth Centered Inertial (ECI) frame. Additionally, change the figure of Earth from sphere to WGS84 ellipsoid, and update the calculations for ground range to use Vincenty's inverse formula. Finally, create the following plots to show the rocket's flight in the ECEF frame:
    * 'Trajectory Metrics (ECEF frame)'
    * 'Trajectory Across Limb of Earth (ECEF frame)'

15) Extend the simulation from 2 dimensions to 3 dimensions, and plot a rocket trajectory on a 3D globe.
16) Implement a launch azimuth by converting from ECI to a topocentric-horizon frame (e.g. SEZ), applying the azimuth rotation to the velocity vector, and then converting it back to ECI. Additionally, rework the pitch program so that the steering vector is applied in the topocentric-horizon frame along the instantaneous azimuth. Finally, show the trajectory of a rocket having a non-zero launch azimuth on a 3D globe and generate a 2D ground trace plot using latitude, longitude, and time data.
17) Compute the six classical orbital elements from the state vector at burnout
18) Simulate the effects of Earth's oblateness on the rocket's motion by implementing J2 gravity
19) Create a GUI for user input and display of outputs, then package the program into an executable using Pyinstaller

## Possible Future Extensions
1) Add ICBM targeting: specify lat and lon of launch point and lat and lon of impact point, and compute ICBM's required altitude, velocity, and flight path angle at final stage burnout
2) Extend simulation from 3-DOF to 6-DOF. Must account for rocket body inertia, CG, CP.
3) Simulate ICBM intercept using interceptor rocket with kill vehicle payload. Use Lambert intercept guidance and assume perfect knowledge of threat state vector at burnout, i.e., the Predicted Intercept Point is the intercept point, so kill vehicle need not maneuver at all to hit the target (Zero Effort Miss distance is zero)
4) Add noise to ICBM state vector, then extend the interceptor's Lambert intercept guidance solution with proportional navigation, i.e., kill vehicle must null the Predicted Intercept Point errors at endgame using its divert and attitude control thrusters


