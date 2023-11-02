#Function to pull NEXRAD Scans from the AWS interface using a netcdf with TC
#supercell cases. 
#The function requires 3 inputs:
#   1. casefile, which is a directory path to the .nc file needed.
#   Ex: casefile = "/Users/colinwelty/Desktop/Python/capstone/capstone_2023.nc"
#   2. case, which is the string name of the desired supercell case.
#   Ex: case = "201708-053"
#   3. radar_id, which is the name of the radar desired for the scan
#   Ex: radar_id = "KHGX"
#   4. track_point, which is the time indice of the mesocyclone
#   Ex: track_point = 0
# 
#Author: Colin Welty (October 2023), with help from nexradaws documentation

from netCDF4 import Dataset, num2date, MFDataset
import numpy as np
from datetime import datetime
import xarray as xr
import tempfile
import nexradaws


def getradarscans(casefile, case, radar_id, track_point):
    conn = nexradaws.NexradAwsInterface()
    templocation = tempfile.mkdtemp()

    # Open the NetCDF file using xarray
    dataset = xr.open_dataset(casefile)

    # Filter the dataset based on the mesocyclone_id
    case_used = dataset.where(dataset.mesocyclone_id == case, drop=True)    

    year = int(case_used['mesocyclone_year'].values[track_point])
    month = int(case_used['mesocyclone_month'].values[track_point])
    day = int(case_used['mesocyclone_day'].values[track_point])
    hours = case_used['mesocyclone_hour'].values
    hour1 = int(hours[track_point])
    hour2 = int(hours[track_point+1])
    hour2+=hour2
    minutes = case_used['mesocyclone_minute'].values
    minute1 = int(minutes[track_point])
    minute2 = int(minutes[track_point+1])
    start = datetime(year,month,day,hour1,minute1)
    end = datetime(year,month,day,hour2,minute2)
    scans = conn.get_avail_scans_in_range(start, end, radar_id)
    results = conn.download(scans, templocation)
    
    # extract mesocyclone center
    lat = case_used['mesocyclone_latitude'].values[track_point]
    lon = case_used['mesocyclone_longitude'].values[track_point]
    
    return results, lat, lon