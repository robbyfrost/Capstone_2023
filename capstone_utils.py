# ------------------------------------------------
# Name: capstone_utils.py
# Author: Robert M. Frost
# University of Oklahoma
# Created: 01 November 2023
# Purpose: Functions to assist with TC supercell 
# capstone.  
# 
# timeheight: Read in radar data from TC supercell 
# case,horizontally average a variable in a column, 
# and combine into timeheight series
# ------------------------------------------------
import pyart
import numpy as np
import xarray as xr
import os
from datetime import datetime
import pandas as pd
import warnings
from PIL import Image
from netCDF4 import Dataset
# ------------------------------------------------
def timeheight(drad, dtrack, meso_id):
    """Read in radar data from TC supercell case, horizontally average a 
    variable in a column, and combine into timeheight series.

    :param str drad: Directory in which nexrad files are stored
    :param str dtrack: Filepath to capstone_2023.nc file
    :param str meso_id: name of case to plot
    """
    # read in tracking dataset
    tracks = xr.open_dataset(dtrack)
    case_meso = tracks.where(tracks.mesocyclone_id == meso_id, drop=True)
    case_meso["mesocyclone_longitude"] = case_meso.mesocyclone_longitude - 360
    # number of track points
    nt = case_meso.mesocyclone_track_points.size
    # Create array to store radar grid
    g_all = np.empty((nt, 241, 241))
    # Create arrays to store grid lat/lon
    glat_all, glon_all = np.empty((nt, 241, 241)), np.empty((nt, 241, 241))
    # Create arrays to store meso location
    mlat_all, mlon_all = np.empty(nt), np.empty(nt)
    # list of radar file names
    rfiles = os.listdir(drad)
    nfile = len(rfiles)
    
    # Loop over files in the directory
    for i in range(nt):
        # read in radar file
        radar = pyart.io.read(f"{drad}{rfiles[i]}")
        # extract radar position
        # radar_lat = radar.latitude['data'][0]
        # radar_lon = radar.longitude['data'][0]
        # mask out last 10 gates of each ray, this removes the "ring" around the radar.
        radar.fields["reflectivity"]["data"][:, -10:] = np.ma.masked
        # exclude masked gates from the gridding
        gatefilter = pyart.filters.GateFilter(radar)
        gatefilter.exclude_transition()
        gatefilter.exclude_masked("reflectivity")
        # grid data to cartesian format
        grid = pyart.map.grid_from_radars(
            (radar,),
            gatefilters=(gatefilter,),
            grid_shape=(1, 241, 241),
            grid_limits=((2000, 2000), (-123000.0, 123000.0), (-123000.0, 123000.0)),
            fields=["reflectivity", "velocity"],
        )
        # add to array of grids
        g_all[i] = grid.fields['reflectivity']['data']
        g_all = g_all[:,np.newaxis,:,:]
        # extract grid lat/lon, put in array
        glat_all[i], glon_all[i] = grid.point_latitude['data'], grid.point_longitude['data']
        # arrays of meso locations
        mlat_all[i], mlon_all[i] = case_meso.mesocyclone_latitude[i], case_meso.mesocyclone_longitude[i]

        # return meso location in x and y position relative to radar
        # meso_x, meso_y = geographic_to_cartesian_aeqd(meso_lon, meso_lat, radar_lon, radar_lat)
    
    # convert lists to arrays
    # mask data around mesocyclone
    mask_grid = mask_tc_track(mlon_all, mlat_all, 1, g_all, glon_all, glat_all)

# -------------------------
def mask_tc_track(clon_tmp, clat, rmax, grid, lon_tmp, lat):
    # Function to apply masking beyond a radius threshold WRT to a tracked
    #   object, e.g., a TC.
    # 
    # INPUTS:
    #       clon_tmp - array of object longitude points
    #       clat - array of object latitude points
    #       rmax - radius limit [degrees lon/lat]
    #       grid - cartesian formatted data
    #       lon_tmp - array of grid longitude points as (x,y) [deg]
    #       lat - " " grid latitude points
    #       t0, t1 - bounding time indices of var, assuming shorter in time than clon/clat
    # 
    # RETURNS:
    #       returns a masked array of identical shape to var.
    # 
    # James Ruppert  
    # jruppert@ou.edu
    # September 2022

    # Input dimensions
    nt,nz,nx1,nx2 = grid[:,0,0,0].size, grid[0,:].size, len(lon_tmp), len(lat)

    clon_tmp = np.ma.masked_invalid(clon_tmp, copy=True)

    # Function to account for crossing of the Intl Date Line
    def dateline_lon_shift(lon_in, reverse):
        if reverse == 0:
            lon_offset = np.ma.zeros(lon_in.shape)
            lon_offset[np.where(lon_in < 0)] += 360
        else:
            lon_offset = np.ma.zeros(lon_in.shape)
            lon_offset[np.where(lon_in > 180)] -= 360
        # return lon_in + lon_offset
        return lon_offset

    # Check for crossing Date Line
    if (lon_tmp.min() < 0) and (lon_tmp.max() > 0):
        lon_offset = dateline_lon_shift(lon_tmp, reverse=0)
        clon_offset = dateline_lon_shift(clon_tmp, reverse=0)
    else:
        lon_offset = 0
        clon_offset = 0
    lon = np.ma.copy(lon_tmp) + lon_offset
    clon = np.ma.copy(clon_tmp) + clon_offset

    # Calculate radius from center as array(time,x,y)
    lon3d = np.repeat(lon[np.newaxis,:,:], nt, axis=0)
    lat3d = np.repeat(lat[np.newaxis,:,:], nt, axis=0)
    lon3d -= clon[:,np.newaxis,np.newaxis]
    lat3d -= clat[:,np.newaxis,np.newaxis]
    radius3d = np.sqrt( lon3d**2 + lat3d**2 )
    radius3d = np.ma.masked_invalid(radius3d, copy=True)

    # Add vertical dimension to match shape of var
    radius4d = np.repeat(radius3d[:,np.newaxis,:,:], nz, axis=1)
    radius4d = np.ma.masked_invalid(radius4d, copy=True)

    # Apply mask
    var_mask = np.ma.masked_where(radius4d > rmax, var, copy=True)

    # Mask out domain edges, within 0.5*r_max from boundaries
    edge = 1 # deg
    npts = int(np.rint(edge/(lon[0,1]-lon[0,0])))
    # var_mask.mask[:,:,0:npts,:] = True
    var_mask.mask[:,:,0:npts,:] = True
    var_mask.mask[:,:,:,0:npts] = True
    var_mask.mask[:,:,-npts:,:] = True
    var_mask.mask[:,:,:,-npts:] = True

    return var_mask

# -------------------------
def geographic_to_cartesian_aeqd(lon, lat, lon_0, lat_0, R=6370997.0):
    """
    FROM PYART Azimuthal equidistant geographic to Cartesian coordinate transform.

    Transform a set of geographic coordinates (lat, lon) to
    Cartesian/Cartographic coordinates (x, y) using a azimuthal equidistant
    map projection [1]_.

    .. math::

        x = R * k * \\cos(lat) * \\sin(lon - lon_0)

        y = R * k * [\\cos(lat_0) * \\sin(lat) -
                     \\sin(lat_0) * \\cos(lat) * \\cos(lon - lon_0)]

        k = c / \\sin(c)

        c = \\arccos(\\sin(lat_0) * \\sin(lat) +
                     \\cos(lat_0) * \\cos(lat) * \\cos(lon - lon_0))

    Where x, y are the Cartesian position from the center of projection;
    lat, lon the corresponding latitude and longitude; lat_0, lon_0 are the
    latitude and longitude of the center of the projection; R is the radius of
    the earth (defaults to ~6371 km).

    Parameters
    ----------
    lon, lat : array-like
        Longitude and latitude coordinates in degrees.
    lon_0, lat_0 : float
        Longitude and latitude, in degrees, of the center of the projection.
    R : float, optional
        Earth radius in the same units as x and y. The default value is in
        units of meters.

    Returns
    -------
    x, y : array
        Cartesian coordinates in the same units as R, typically meters.

    References
    ----------
    .. [1] Snyder, J. P. Map Projections--A Working Manual. U. S. Geological
        Survey Professional Paper 1395, 1987, pp. 191-202.

    """
    lon = np.atleast_1d(np.asarray(lon))
    lat = np.atleast_1d(np.asarray(lat))

    lon_rad = np.deg2rad(lon)
    lat_rad = np.deg2rad(lat)

    lat_0_rad = np.deg2rad(lat_0)
    lon_0_rad = np.deg2rad(lon_0)

    lon_diff_rad = lon_rad - lon_0_rad

    # calculate the arccos after ensuring all values in valid domain, [-1, 1]
    arg_arccos = np.sin(lat_0_rad) * np.sin(lat_rad) + np.cos(lat_0_rad) * np.cos(
        lat_rad
    ) * np.cos(lon_diff_rad)
    arg_arccos[arg_arccos > 1] = 1
    arg_arccos[arg_arccos < -1] = -1
    c = np.arccos(arg_arccos)
    with warnings.catch_warnings():
        # division by zero may occur here but is properly addressed below so
        # the warnings can be ignored
        warnings.simplefilter("ignore", RuntimeWarning)
        k = c / np.sin(c)
    # fix cases where k is undefined (c is zero), k should be 1
    k[c == 0] = 1

    x = R * k * np.cos(lat_rad) * np.sin(lon_diff_rad)
    y = (
        R
        * k
        * (
            np.cos(lat_0_rad) * np.sin(lat_rad)
            - np.sin(lat_0_rad) * np.cos(lat_rad) * np.cos(lon_diff_rad)
        )
    )

    return x, y

# -------------------------
def column_vertical_profile(
    radar, latitude, longitude, azimuth_spread=3, spatial_spread=3
):
    """
    Given the location (in latitude, longitude) of a target, return the rays
    that correspond to radar column above the target, allowing for user
    defined range of azimuths and range gates to be included within this
    extraction.

    Parameters
    ----------
    radar : pyart.core.Radar Object
        Py-ART Radar Object from which distance to the target, along
        with gates above the target, will be calculated.
    latitude : float, [degrees]
        Latitude, in degrees North, of the target.
    longitude : float, [degrees]
        Longitude, in degrees East, of the target.
    azimuth_spread : int
        Number of azimuth angles to include within extraction list
    spatial_range : int
        Number of range gates to include within the extraction

    Function Calls
    --------------
    sphere_distance
    for_azimuth
    get_sweep_rays
    subset_fields
    assemble_column

    Returns
    -------
    column : xarray
        Xarray Dataset containing the radar column above the target for
        the various fields within the radar object.

    References
    ----------
    Murphy, A. M., A. Ryzhkov, and P. Zhang, 2020: Columnar Vertical
    Profile (CVP) Methodology for Validating Polarimetric Radar Retrievals
    in Ice Using In Situ Aircraft Measurements. J. Atmos. Oceanic Technol.,
    37, 1623–1642, https://doi.org/10.1175/JTECH-D-20-0011.1.

    Bukovčić, P., A. Ryzhkov, and D. Zrnić, 2020: Polarimetric Relations for
    Snow Estimation—Radar Verification. J. Appl. Meteor. Climatol.,
    59, 991–1009, https://doi.org/10.1175/JAMC-D-19-0140.1.
    """

    # Define the spatial range to use within extraction
    spatial_range = radar.range["meters_between_gates"] * spatial_spread

    # Define a dictionary structure to contain the extracted features
    total_moment = {key: [] for key in radar.fields.keys()}
    total_moment.update({"height": [], "time_offset": []})

    # Define the start of the radar volume
    base_time = pd.to_datetime(radar.time["units"][14:]).to_numpy()

    # call the sphere_distance function
    dis = sphere_distance(
        radar.latitude["data"][0], latitude, radar.longitude["data"][0], longitude
    )
    # calculate forward azimuth angle
    forazi = for_azimuth(
        radar.latitude["data"][0], latitude, radar.longitude["data"][0], longitude
    )

    # Iterate through radar sweeps, extract desired section
    for sweep in radar.iter_slice():
        moment = {key: [] for key in radar.fields.keys()}
        zgates = []
        gate_time = []

        # call the new sweep rays
        center, spread = get_sweep_rays(
            radar.azimuth["data"][sweep], forazi, azimuth_spread=azimuth_spread
        )

        # add the start indice of each ray
        center = [x + sweep.start for x in center]
        spread = [x + sweep.start for x in spread]
        # Correct the spread indices to remove the centerline
        spread = [x for x in spread if x not in center]

        # For the ray(s) directly over the target, extract and average fields
        for ray in center:
            # Convert gates from antenna or cartesian coordinates
            (rhi_x, rhi_y, rhi_z) = antenna_vectors_to_cartesian(
                radar.range["data"],
                radar.azimuth["data"][ray],
                radar.elevation["data"][ray],
                edges=False,
            )
            # Calculate distance to target
            rhidis = np.sqrt((rhi_x**2) + (rhi_y**2)) * np.sign(rhi_z)
            # Calculate target gate
            tar_gate = np.nonzero(np.abs(rhidis[0, :] - dis) < spatial_range)[
                0
            ].tolist()
            # Subset the radar fields for the target locations
            subset = subset_fields(radar, ray, tar_gate)
            # Add back to the total dictionary
            moment = {key: moment[key] + subset[key] for key in moment}
            # Add radar elevation to height gates
            # to define height as center of each gate above sea level
            zgates.append(np.ma.mean(rhi_z[0, tar_gate] + radar.altitude["data"][0]))
            # Determine the time for the individual gates
            gate_time.append(radar.time["data"][ray])

        # Convert to Cartesian Coordinates
        # Determine the center of each gate for the subsetted rays.
        for ray in spread:
            (rhi_x, rhi_y, rhi_z) = antenna_vectors_to_cartesian(
                radar.range["data"],
                radar.azimuth["data"][ray],
                radar.elevation["data"][ray],
                edges=False,
            )
            # Calculate distance to target
            rhidis = np.sqrt((rhi_x**2) + (rhi_y**2)) * np.sign(rhi_z)
            # Calculate target gate
            tar_gate = np.nonzero(np.abs(rhidis[0, :] - dis) < spatial_range)[
                0
            ].tolist()
            # Subset the radar fields for the target locations
            subset = subset_fields(radar, ray, tar_gate)
            # Add back to the sweep dictionary
            moment = {key: moment[key] + subset[key] for key in moment}
            # Add radar elevation to height gates
            # to define height as center of each gate above sea level
            zgates.append(np.ma.mean(rhi_z[0, tar_gate] + radar.altitude["data"][0]))
            # Determine the time for the individual gates
            gate_time.append(radar.time["data"][ray])

        # Average all azimuth moments into a single value for the sweep
        for key in total_moment:
            if key == "height":
                total_moment[key].append(np.ma.mean(np.ma.masked_invalid(zgates)))
            elif key == "time_offset":
                total_moment[key].append(np.round(np.ma.mean(np.array(gate_time)), 4))
            else:
                total_moment[key].append(
                    np.round(np.ma.mean(np.ma.masked_invalid(moment[key])), 4)
                )

    # Add the base time for the radar
    total_moment.update({"base_time": base_time})

    # Convert to xarray
    return assemble_column(radar, total_moment, forazi, dis, latitude, longitude)

# -------------------------
def make_gif(dfig, dout):
    """Append figures into a gif and output

    :param str dfig: Directory in which figures lie
    :param str dout: File path and name for gif
    """
    # create list to store figure paths
    image_all = []
    
    # list of file names in dfig
    files = os.listdir(dfig)
    # sorted alphabetically
    sorted_files = sorted(files)
    
    # loop over files
    for i in range(len(files)):
        image_path = f"{dfig}{sorted_files[i]}"
        image_all.append(image_path)

    # Create a list to store the image objects
    images = []

    # Open and append each image to the list
    for file in image_all:
        image = Image.open(file)
        images.append(image)

    # Save the images as an animated GIF
    images[0].save(dout, save_all=True, append_images=images[1:], optimize=False, duration=1000, loop=0)

    print(f"GIF saved successfully at {dout}")

# -------------------------
# test function
drad = "/Users/robertfrost/Documents/boomsoon/2023_semester_2/capstone/analysis/Capstone_2023/data/lightning_case1/"
dtrack = "/Users/robertfrost/Documents/boomsoon/2023_semester_2/capstone/analysis/capstone_2023.nc"
meso_id = "201708-053"
timeheight(drad, dtrack, meso_id)