#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 16:43:58 2023

@author: robertfrost
"""

import numpy as np
import xarray as xr
import pyart
from getradarscans import getradarscans
import matplotlib.pyplot as plt
import datetime
import glob

casefile_path = "/Users/robertfrost/Documents/boomsoon/2023_semester_2/capstone/data/capstone_2023.nc"
case_str = "201708-053"
radar_id = "KHGX"
track_point = 0

casefile = xr.open_dataset(casefile_path)
case = casefile.where(casefile.mesocyclone_id == case_str, drop=True)

year = int(case['mesocyclone_year'].values[track_point])
month = int(case['mesocyclone_month'].values[track_point])
if month < 10:
    month = f"0{month}"
day = int(case['mesocyclone_day'].values[track_point])
if day < 10:
    day = f"0{day}"
hour = int(case['mesocyclone_hour'].values[track_point])
if hour < 10:
    hour = f"0{hour}"
minute = int(case['mesocyclone_minute'].values[track_point])
if minute < 10:
    minute = f"0{minute}"

radar_path = f"{casefile_path[:68]}lightning_case1/{radar_id}{year}{month}{day}_{minute}[0-60]*"

radar = pyart.io.read(glob.glob(radar_path[0]))

# casefile = "/Users/robertfrost/Documents/boomsoon/2023_semester_2/capstone/data/capstone_2023.nc"
# case = "201708-053"
# radar_id = "KHGX"
# track_point = 0

# results, lat, lon = getradarscans(casefile, case, radar_id, track_point)

# for i, scan in enumerate(results.iter_success()):
#     radar = scan.open_pyart()
#     display = pyart.graph.RadarMapDisplay(radar)
    
#     center_lat = lat
#     center_lon = lon
    
#     # Specify the radius in degrees (adjust as needed)
#     radius_degrees = 0.10
#     # Calculate the lat/lon boundaries based on the center and radius
#     lat_min = center_lat - radius_degrees
#     lat_max = center_lat + radius_degrees
#     lon_min = center_lon - radius_degrees
#     lon_max = center_lon + radius_degrees
    
#     # fig, ax = plt.subplots(figsize=(12,10))
#     display.plot_ppi_map("velocity", sweep=4, vmin=-20, vmax=20,
#                           min_lon=lon_min, max_lon=lon_max, 
#                           min_lat=lat_min, max_lat=lat_max)
    
#     # ax.set_aspect('equal')
    
#     # plt.show()
#     # plt.close()