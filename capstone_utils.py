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
# ------------------------------------------------
# def timeheight(drad, dtrack, meso_id):
#     """Read in radar data from TC supercell case, horizontally average a 
#     variable in a column, and combine into timeheight series.

#     :param str drad: Directory in which nexrad files are stored
#     :param str dtrack: Filepath to capstone_2023.nc file
#     :param str meso_id: name of case to plot
#     """
    
#     # Create lists to store wspd, u, and v DataArrays
#     refl_all = []
#     # Create a list to store the timestamps
#     times = []
    
#     # Loop over files in the directory
#     for i in range(os.listdir(drad)):
#         # read in radar file
#         radar = pyart.io.read(f"{drad}{filename}")
#         # read in tracking dataset
#         tracks = xr.open_dataset(dtrack)
#         case_meso = xr.where(tracks.mesocyclone_id == meso_id, drop=True)

#         # meso lat/lon
#         lat = case_meso.mesocyclone_latitude[]

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