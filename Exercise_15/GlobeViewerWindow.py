# -*- coding: utf-8 -*-
import sys
import pyvista as pv
from pyvistaqt import QtInteractor
from PySide6.QtWidgets import (
    QApplication, QMainWindow, QWidget, 
    QVBoxLayout, QHBoxLayout, QSlider, 
    QPushButton, QStyle, QLabel, QFrame
)
from PySide6.QtCore import Qt, QTimer
import numpy as np
import geophysical_tools as gt
from symbols import aries_symbol
import pickle
from symbols import deg_symbol
import matplotlib.colors as mcolors

font_file_path = r".\Misc\DejaVuSans.ttf"
globe_map_file_path = r".\Misc\blue_marble_3.jpg"
mm3_rocket_file_path=r".\Misc\mm3_rocket.pkl"
f9_rocket_file_path=r".\Misc\f9_rocket.pkl"
v2_rocket_file_path=r".\Misc\v2_rocket.pkl"

class GlobeViewerWindow(QMainWindow):
    def __init__(self, Rocket):
        super().__init__()
        self.setWindowTitle("Globe Viewer Window")
        self.Rocket = Rocket
        
        # Central Widget & Outer Horizontal Layout (Spans full height)
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_h_layout = QHBoxLayout(central_widget)
        main_h_layout.setContentsMargins(0, 0, 0, 0)  # Optional: remove exterior gaps

        # 2. Sidebar Widget (Extends full vertical height on the left)
        self.sidebar_widget = QWidget()
        self.sidebar_layout = QVBoxLayout(self.sidebar_widget) # Changed to Vertical for sidebar items
        self.sidebar_widget.setFixedWidth(340)
        
        # --- sidebar controls ---
        title_label = QLabel(f"{self.Rocket.config['name']} Rocket Trajectory Simulation")
        title_label.setStyleSheet("font-size: 28px; font-weight: bold;")
        title_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        title_label.setWordWrap(True)
        self.sidebar_layout.addWidget(title_label)
        self.sidebar_layout.addSpacing(5)
        
        # Create a line frame
        line = QFrame()
        line.setFrameShape(QFrame.Shape.HLine)  # Use VLine for vertical side-by-side layouts
        line.setFrameShadow(QFrame.Shadow.Sunken)
        self.sidebar_layout.addWidget(line)
        
        self.sidebar_layout.addSpacing(5)
        time_info_label = QLabel("Time Info:")
        time_info_label.setStyleSheet("font-size: 24px; font-weight: bold;")
        self.sidebar_layout.addWidget(time_info_label)
        
        ut1datetime_item_layout = QHBoxLayout()
        ut1datetime_item_layout.setContentsMargins(15, 0, 15, 0)
        ut1datetime_text_label = QLabel("UT1:")
        ut1datetime_text_label.setStyleSheet("font-size: 20px")
        self.ut1datetime_value_label = QLabel()
        self.ut1datetime_value_label.setStyleSheet("font-size: 20px")
        ut1datetime_item_layout.addWidget(ut1datetime_text_label)
        ut1datetime_item_layout.addStretch()
        ut1datetime_item_layout.addWidget(self.ut1datetime_value_label)
        self.sidebar_layout.addLayout(ut1datetime_item_layout)
        
        talo_item_layout = QHBoxLayout()
        talo_item_layout.setContentsMargins(15, 0, 15, 0)
        talo_text_label = QLabel("TALO:")
        talo_text_label.setStyleSheet("font-size: 20px")
        self.talo_value_label = QLabel()
        self.talo_value_label.setStyleSheet("font-size: 20px")
        talo_item_layout.addWidget(talo_text_label)
        talo_item_layout.addStretch()
        talo_item_layout.addWidget(self.talo_value_label)
        self.sidebar_layout.addLayout(talo_item_layout)
        
        self.sidebar_layout.addSpacing(10)
        self.kinematic_info_label = QLabel("Kinematics:")
        self.kinematic_info_label.setStyleSheet("font-size: 24px; font-weight: bold;")
        self.sidebar_layout.addWidget(self.kinematic_info_label)
        
        lat_item_layout = QHBoxLayout()
        lat_item_layout.setContentsMargins(15, 0, 15, 0)
        lat_text_label = QLabel("Latitude:")
        lat_text_label.setStyleSheet("font-size: 20px;")
        self.lat_value_label = QLabel()
        self.lat_value_label.setStyleSheet("font-size: 20px;")
        lat_item_layout.addWidget(lat_text_label)
        lat_item_layout.addStretch()
        lat_item_layout.addWidget(self.lat_value_label)
        self.sidebar_layout.addLayout(lat_item_layout)
        
        lon_item_layout = QHBoxLayout()
        lon_item_layout.setContentsMargins(15, 0, 15, 0)
        lon_text_label = QLabel("Longitude:")
        lon_text_label.setStyleSheet("font-size: 20px;")
        self.lon_value_label = QLabel()
        self.lon_value_label.setStyleSheet("font-size: 20px;")
        lon_item_layout.addWidget(lon_text_label)
        lon_item_layout.addStretch()
        lon_item_layout.addWidget(self.lon_value_label)
        self.sidebar_layout.addLayout(lon_item_layout)
        
        alt_item_layout = QHBoxLayout()
        alt_item_layout.setContentsMargins(15, 0, 15, 0)
        alt_text_label = QLabel("Altitude:")
        alt_text_label.setStyleSheet("font-size: 20px;")
        self.alt_value_label = QLabel()
        self.alt_value_label.setStyleSheet("font-size: 20px;")
        alt_item_layout.addWidget(alt_text_label)
        alt_item_layout.addStretch()
        alt_item_layout.addWidget(self.alt_value_label)
        self.sidebar_layout.addLayout(alt_item_layout)
        
        ground_range_item_layout = QHBoxLayout()
        ground_range_item_layout.setContentsMargins(15, 0, 15, 0)
        ground_range_text_label = QLabel("Ground range:")
        ground_range_text_label.setStyleSheet("font-size: 20px;")
        self.ground_range_value_label = QLabel()
        self.ground_range_value_label.setStyleSheet("font-size: 20px;")
        ground_range_item_layout.addWidget(ground_range_text_label)
        ground_range_item_layout.addStretch()
        ground_range_item_layout.addWidget(self.ground_range_value_label)
        self.sidebar_layout.addLayout(ground_range_item_layout)
               
        ecef_total_vel_item_layout = QHBoxLayout()
        ecef_total_vel_item_layout.setContentsMargins(15, 0, 15, 0)
        ecef_total_vel_text_label = QLabel("v<sub>total,ECEF</sub>: ")
        ecef_total_vel_text_label.setStyleSheet("font-size: 20px;")
        self.ecef_total_vel_value_label = QLabel()
        self.ecef_total_vel_value_label.setStyleSheet("font-size: 20px;")
        ecef_total_vel_item_layout.addWidget(ecef_total_vel_text_label)
        ecef_total_vel_item_layout.addStretch()
        ecef_total_vel_item_layout.addWidget(self.ecef_total_vel_value_label)
        self.sidebar_layout.addLayout(ecef_total_vel_item_layout)
        
        eci_total_vel_item_layout = QHBoxLayout()
        eci_total_vel_item_layout.setContentsMargins(15, 0, 15, 0)
        eci_total_vel_text_label = QLabel("v<sub>total,ECI</sub>: ")
        eci_total_vel_text_label.setStyleSheet("font-size: 20px;")
        self.eci_total_vel_value_label = QLabel()
        self.eci_total_vel_value_label.setStyleSheet("font-size: 20px;")
        eci_total_vel_item_layout.addWidget(eci_total_vel_text_label)
        eci_total_vel_item_layout.addStretch()
        eci_total_vel_item_layout.addWidget(self.eci_total_vel_value_label)
        self.sidebar_layout.addLayout(eci_total_vel_item_layout)
        
        self.sidebar_layout.addSpacing(10)
        flight_phases_title_label = QLabel("Flight Phases:")
        flight_phases_title_label.setStyleSheet("font-size: 24px; font-weight: bold;")
        self.sidebar_layout.addWidget(flight_phases_title_label)
        
        # Build dynamic legend items for each flight phase
        self.phase_labels = []
        for i, phase_name in enumerate(self.Rocket.phase_t_idx_dict.keys()):
            # Get color from Rocket phase_color_dict
            color = self.Rocket.phase_color_dict[i]
            rgba = mcolors.to_rgba(color)
            r,g,b = [int(x*255) for x in rgba[:3]]
            css_color = f"rgb({r}, {g}, {b})"

            # Horizontal layout for the colored box + phase text
            legend_item_layout = QHBoxLayout()
            legend_item_layout.setContentsMargins(15, 0, 15, 0)

            # Colored color-swatch icon
            color_box = QFrame()
            color_box.setFixedSize(18, 18)
            color_box.setStyleSheet(f"background-color: {css_color}; border-radius: 3px; border: 1px solid #333;")

            # Phase label
            phase_label = QLabel(phase_name)
            phase_label.setStyleSheet("font-size: 20px; font-weight: normal; color: gray;")
            phase_label.setWordWrap(True)
            self.phase_labels.append(phase_label)

            legend_item_layout.addWidget(color_box)
            legend_item_layout.addWidget(phase_label)
            legend_item_layout.addStretch()  # Keep left-aligned

            self.sidebar_layout.addLayout(legend_item_layout)
        
        # push all added widgets to the top
        self.sidebar_layout.addStretch()
        
        # Main Content Container (Right side: Viewport + Bottom Bar stacked vertically)
        content_widget = QWidget()
        content_v_layout = QVBoxLayout(content_widget)
        content_v_layout.setContentsMargins(0, 0, 0, 0)

        # Bottom bar widget (Horizontal layout for Play/pause button, FF button, and Slider)
        bottombar_layout = QVBoxLayout()
        bottombar_widget = QWidget()
        bottombar_widget.setLayout(bottombar_layout)
        bottombar_widget.setFixedHeight(50)

        # Embed PyVista Viewport into the right content layout
        self.plotter = QtInteractor(content_widget)
        content_v_layout.addWidget(self.plotter.interactor)
        content_v_layout.addWidget(bottombar_widget)

        # Add Sidebar and Content Area to the Root Horizontal Layout
        main_h_layout.addWidget(self.sidebar_widget)
        main_h_layout.addWidget(content_widget)

        # Controls Container 
        controls_layout = QHBoxLayout()
        
        # --- Play / Pause Button Setup ---
        self.play_icon = self.style().standardIcon(QStyle.StandardPixmap.SP_MediaPlay)
        self.pause_icon = self.style().standardIcon(QStyle.StandardPixmap.SP_MediaPause)

        self.play_button = QPushButton()
        self.play_button.setIcon(self.play_icon)
        self.play_button.setFixedWidth(40)
        self.play_button.clicked.connect(self.toggle_play)
        controls_layout.addWidget(self.play_button)
        
        # --- Speed & Fast-Forward Setup ---
        self.speeds = [1, 5, 10, 50, 100]
        self.speed_idx = 0  # 0 corresponds to 1x speed

        # --- Fast-Forward Button ---
        self.ff_button = QPushButton("1x")
        self.ff_button.setFixedWidth(50)
        self.ff_button.setStyleSheet("font-size: 14px;")
        self.ff_button.clicked.connect(self.cycle_speed)
        controls_layout.addWidget(self.ff_button)

        # Constant label properties
        self.axes_label_distance = gt.R_eq * 1.6
        self.axes_label_size = 20
                
        # Build Scene
        self.render_scene(Rocket)
                
        # --- Time Slider Setup ---
        self.initial_slider_time_idx = 0
        self.time_slider = QSlider(Qt.Orientation.Horizontal)
        self.time_slider.setMinimum(0)
        self.time_slider.setMaximum(len(Rocket.times) - 1)
        self.time_slider.setValue(self.initial_slider_time_idx)
        controls_layout.addWidget(self.time_slider)
        
        bottombar_layout.addLayout(controls_layout)

        # --- Animation Timer Setup ---
        self.timer = QTimer(self)
        self.timer.setInterval(33)  # Fixed ~30 FPS delay
        self.timer.timeout.connect(self.advance_slider)

        #initialize the rotation of objects in the scene and the time label
        self.time_slider.valueChanged.connect(self.on_time_idx_change)
        
        # Apply rotation of ECEF elements by the GMST corresponding to the initial time index
        self.on_time_idx_change(self.initial_slider_time_idx)

    # --- Animation Control Methods ---
    def cycle_speed(self):
        # Cycle through index: 0 (1x) -> 1 (2x) -> 2 (3x) -> 0 (1x)
        self.speed_idx = (self.speed_idx + 1) % len(self.speeds)
        current_speed = self.speeds[self.speed_idx]
        
        # Update button text
        self.ff_button.setText(f"{current_speed}x")
        
    def toggle_play(self):
            if self.timer.isActive():
                self.timer.stop()
                self.play_button.setIcon(self.play_icon)
            else:
                if self.time_slider.value() >= self.time_slider.maximum():
                    self.time_slider.setValue(0)
                
                self.timer.start()
                self.play_button.setIcon(self.pause_icon)

    def advance_slider(self):
        current_val = self.time_slider.value()
        step = self.speeds[self.speed_idx]
        
        if current_val + step <= self.time_slider.maximum():
            self.time_slider.setValue(current_val + step)
        else:
            # Snap to final frame and pause when reaching the end
            self.time_slider.setValue(self.time_slider.maximum())
            self.timer.stop()
            self.play_button.setIcon(self.play_icon)

    def on_time_idx_change(self, time_idx: int):
        # --- Dynamic Trajectory Update ---  
        for i, key in enumerate(self.Rocket.phase_t_idx_dict.keys()):
            start_idx = self.Rocket.phase_t_idx_dict[key][0]
            stop_idx = self.Rocket.phase_t_idx_dict[key][1]
            polydata = self.traj_polydata_list[i]
            actor = self.actor_traj_phases[i]
    
            # --- Dynamic Phase Bolding ---
            is_active_phase = start_idx <= time_idx <= stop_idx
            if is_active_phase:
                self.phase_labels[i].setStyleSheet("font-size: 20px; font-weight: bold; color: black;")
            else:
                self.phase_labels[i].setStyleSheet("font-size: 20px; font-weight: normal; color: gray;")

            # --- Trajectory Segment Visibility & Updates ---
            
            if time_idx < start_idx:
                # Phase hasn't started yet -> Hide segment
                actor.SetVisibility(False)
            else:
                # Phase has started or completed -> Show segment
                actor.SetVisibility(True)
                
                # Clamp the current trajectory slice to the end of this phase
                current_stop = min(time_idx, stop_idx)
                
                # Slice active trajectory points up to current time
                active_points = self.Rocket.r_ecef[start_idx : current_stop + 1]
                
                # Handle edge case where phase just started (need at least 2 points for a line)
                if len(active_points) == 1:
                    active_points = np.vstack([active_points, active_points])
                    
                # Re-generate polydata line geometry and copy into the existing mesh object
                updated_line = pv.lines_from_points(active_points)
                polydata.deep_copy(updated_line)
        
        # Update Sidebar Text
        self.ut1datetime_value_label.setText(f'{self.Rocket.ut1datetimes[time_idx]}')
        self.talo_value_label.setText(f't+{self.Rocket.times[time_idx]}[s]')
        self.lat_value_label.setText(f'{round(self.Rocket.lla[time_idx][0],2)}{deg_symbol}')
        self.lon_value_label.setText(f'{round(self.Rocket.lla[time_idx][1],2)}{deg_symbol}')
        self.alt_value_label.setText(f'{round(self.Rocket.lla[time_idx][2]/1000,2)}[km]')
        self.ground_range_value_label.setText(f'{round(self.Rocket.ground_ranges_km[time_idx],2)}[km]')
        self.ecef_total_vel_value_label.setText(f"{round(np.sqrt(self.Rocket.v_ecef[time_idx][0]**2 + self.Rocket.v_ecef[time_idx][1]**2 + self.Rocket.v_ecef[time_idx][2]**2),2)}[m/s]")
        self.eci_total_vel_value_label.setText(fr"{round(np.sqrt(self.Rocket.v[time_idx][0]**2 + self.Rocket.v[time_idx][1]**2 + self.Rocket.v[time_idx][2]**2),2)}[m/s]")

        angle_rad = self.Rocket.GMSTs[time_idx]
        cos_a, sin_a = np.cos(angle_rad), np.sin(angle_rad)

        # 4x4 Homogeneous Transformation Matrix for GPU UserMatrix
        R_z_4x4 = np.array([
            [cos_a, -sin_a, 0, 0],
            [sin_a,  cos_a, 0, 0],
            [0,      0,     1, 0],
            [0,      0,     0, 1]
        ], dtype=np.float64)

        # 2. ZERO-COPY GPU ROTATION for 3D Geometry
        self.actor_spheroid.user_matrix = R_z_4x4

        for actor in self.actor_longitude_lines:
            actor.user_matrix = R_z_4x4

        self.actor_x_arrow.user_matrix = R_z_4x4
        self.actor_y_arrow.user_matrix = R_z_4x4

        for actor in self.actor_traj_phases:
            actor.user_matrix = R_z_4x4

        # 3. Calculate rotated 3D coordinates for X and Y labels
        R_z_3x3 = np.array([
            [cos_a, -sin_a, 0],
            [sin_a,  cos_a, 0],
            [0,      0,     1]
        ])

        new_x_pt = R_z_3x3 @ np.array([self.axes_label_distance, 0, 0])
        new_y_pt = R_z_3x3 @ np.array([0, self.axes_label_distance, 0])

        # Overwrite label overlays
        self.plotter.add_point_labels(
            points=new_x_pt, 
            labels=[r'$X_{{ECEF}}$'],
            name='x_ecef_label',
            font_file=font_file_path,
            font_size=self.axes_label_size, 
            show_points=False,
            shape_opacity=0,
            justification_horizontal='center',
            justification_vertical='center'
        )

        self.plotter.add_point_labels(
            points=new_y_pt, 
            labels=[r'$Y_{{ECEF}}$'],
            name='y_ecef_label',
            font_file=font_file_path,
            font_size=self.axes_label_size, 
            show_points=False,
            shape_opacity=0,
            justification_horizontal='center',
            justification_vertical='center'
        )

        self.plotter.update()
       
    def render_scene(self, Rocket):
        
        # --- Spheroid Mesh ---
        lon_res, lat_res = 24, 24
        lon = np.linspace(0, 2 * np.pi, lon_res)
        lat = np.linspace(-np.pi/2, np.pi/2, lat_res)
        lon_grid, lat_grid = np.meshgrid(lon, lat)

        x = gt.R_eq * np.cos(lat_grid) * np.cos(lon_grid)
        y = gt.R_eq * np.cos(lat_grid) * np.sin(lon_grid)
        z = gt.R_p * np.sin(lat_grid)

        grid = pv.StructuredGrid(x, y, z)
        spheroid = grid.extract_surface(algorithm='dataset_surface')

        u = (lon_grid.flatten() + np.pi) / (2 * np.pi)
        v = (lat_grid.flatten() + np.pi/2) / np.pi

        spheroid.active_texture_coordinates = np.column_stack((u, v))
        texture = pv.read_texture(globe_map_file_path)

        self.plotter.set_background("white")
        self.actor_spheroid = self.plotter.add_mesh(
            spheroid, 
            texture=texture, 
            smooth_shading=True, 
            name="globe_render",
            ambient=0.3,
            diffuse=1,
            specular=0
        )

        # --- Lines of Longitude ---
        lon_spacing_deg = 30
        longitudes = np.arange(0, np.pi*2, lon_spacing_deg*np.pi/180)
        lat_steps = np.linspace(-np.pi/2, np.pi/2, 24)
        self.actor_longitude_lines = []
        for lon_val in longitudes:
            lon_x = gt.R_eq * np.cos(lat_steps) * np.cos(lon_val)
            lon_y = gt.R_eq * np.cos(lat_steps) * np.sin(lon_val)
            lon_z = gt.R_p * np.sin(lat_steps)
            longitude_line_points = np.column_stack((lon_x, lon_y, lon_z))
            longitude_line = pv.MultipleLines(points=longitude_line_points)
            actor = self.plotter.add_mesh(longitude_line, color="gray", line_width=1)
            self.actor_longitude_lines.append(actor)

        # --- Lines of Latitude ---
        lat_spacing_deg = 30
        latitudes = np.arange(0, np.pi*2, lat_spacing_deg*np.pi/180)
        lon_steps = np.linspace(-np.pi/2, np.pi/2, 24)
        for lat_val in latitudes:
            lat_x = gt.R_eq * np.cos(lat_val) * np.cos(lon_steps)
            lat_y = gt.R_eq * np.cos(lat_val) * np.sin(lon_steps)
            lat_z = gt.R_p * np.sin(lat_val) * np.ones(len(lon_steps))
            latitude_line_points = np.column_stack((lat_x, lat_y, lat_z))
            latitude_line = pv.MultipleLines(points=latitude_line_points)
            self.plotter.add_mesh(latitude_line, color="gray", line_width=1)
            
        # --- ECEF Axes ---
        axes_arrow_scale = gt.R_eq * 1.5
        axes_arrow_shaft_radius = 0.005
        axes_arrow_tip_radius = 0.015
        axes_arrow_tip_length = 0.04
        
        x_ecef_arrow = pv.Arrow(
            start=(0,0,0), direction=(1,0,0), scale=axes_arrow_scale,
            shaft_radius=axes_arrow_shaft_radius, tip_radius=axes_arrow_tip_radius, tip_length=axes_arrow_tip_length
        )
        self.actor_x_arrow = self.plotter.add_mesh(x_ecef_arrow, color="red")
        
        y_ecef_arrow = pv.Arrow(
            start=(0,0,0), direction=(0,1,0), scale=axes_arrow_scale,
            shaft_radius=axes_arrow_shaft_radius, tip_radius=axes_arrow_tip_radius, tip_length=axes_arrow_tip_length
        )
        self.actor_y_arrow = self.plotter.add_mesh(y_ecef_arrow, color="green")
        
        z_ecef_arrow = pv.Arrow(
            start=(0,0,0), direction=(0,0,1), scale=axes_arrow_scale,
            shaft_radius=axes_arrow_shaft_radius, tip_radius=axes_arrow_tip_radius, tip_length=axes_arrow_tip_length
        )
        self.plotter.add_mesh(z_ecef_arrow, color="blue")

        self.plotter.add_point_labels(
            points=[0, 0, self.axes_label_distance], 
            labels=[r'$Z_{{ECEF}}$'],
            name='z_ecef_label',
            font_file=font_file_path,
            font_size=self.axes_label_size, 
            show_points=False,
            shape_opacity=0,
            justification_horizontal='center',
            justification_vertical='center'
        )
            
        # --- Vernal Equinox Vector ---
        vernal_eq_arrow = pv.Arrow(
            start=(0,0,0), direction=(1,0,0), scale=axes_arrow_scale,
            shaft_radius=axes_arrow_shaft_radius, tip_radius=axes_arrow_tip_radius, tip_length=axes_arrow_tip_length
        )
        self.plotter.add_mesh(vernal_eq_arrow, color="magenta")
        self.plotter.add_point_labels(
            points=[self.axes_label_distance, 0, 0], 
            labels=[aries_symbol],
            font_file=font_file_path,
            font_size=self.axes_label_size, 
            show_points=False,
            shape_opacity=0,
            justification_horizontal='center',
            justification_vertical='center'
        )
        
        # --- Trajectory Segments Setup ---
        self.actor_traj_phases = []
        self.traj_polydata_list = [] # Store VTK PolyData objects to update points live
    
        for i, key in enumerate(Rocket.phase_t_idx_dict.keys()):
            start_idx = Rocket.phase_t_idx_dict[key][0]
            # Initialize each segment with its first point to prevent rendering errors
            init_pts = Rocket.r_ecef[start_idx:start_idx+1]
            
            # Create line polydata
            polydata = pv.lines_from_points(init_pts)
            
            # Add mesh to scene with render_lines_as_tubes=True
            actor = self.plotter.add_mesh(
                polydata,
                render_lines_as_tubes=True,
                line_width=5,
                color=self.Rocket.phase_color_dict[i],
            )
        
            self.actor_traj_phases.append(actor)
            self.traj_polydata_list.append(polydata)
        
        # Setup Initial View
        # Calculate distance to fit scene
        distance = gt.R_eq * 7.0 
        
        # Position along X looking at origin, with Z as UP vector
        self.plotter.camera_position = [
            (distance, 0, 0),  # Camera Position
            (0, 0, 0),         # Focal Point
            (0, 0, 1)          # View Up Vector (+Z is North)
        ]
        
        #Set camera azimuth to the trajectory's midpoint longitude
        midpoint_lon_deg = (self.Rocket.lla[0,1] + self.Rocket.lla[-1,1]) / 2
        camera_az_deg = midpoint_lon_deg + np.degrees(self.Rocket.GMSTs[0])
        
        #Set camera elevation to the trajectory's midpoint latitude
        camera_el_deg = (self.Rocket.lla[0,0] + self.Rocket.lla[-1,0]) / 2
        
        #Apply camera azimuth and elevation
        self.plotter.camera.azimuth = camera_az_deg
        self.plotter.camera.elevation = camera_el_deg
        self.plotter.update()

    def closeEvent(self, event):
        self.timer.stop()
        self.plotter.close()
        super().closeEvent(event)

if __name__ == '__main__':
    with open(v2_rocket_file_path, "rb") as file:
        v2 = pickle.load(file)
        
    with open(f9_rocket_file_path, "rb") as file:
        f9 = pickle.load(file)
        
    with open(mm3_rocket_file_path, "rb") as file:
        mm3 = pickle.load(file)
        
    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
        
    window = GlobeViewerWindow(mm3)
    
    screen = app.primaryScreen()
    if screen:
        window.setGeometry(screen.availableGeometry())
    window.showMaximized()
    
    sys.exit(app.exec()) if not sys.flags.interactive else app.exec()