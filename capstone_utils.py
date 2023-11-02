#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  1 15:42:41 2023

@author: robertfrost
"""

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
# ------------------------------------------------
def timeheight(drad):
    """Read in radar data from TC supercell case,horizontally average a 
    variable in a column, and combine into timeheight series.

    :param str drad: Directory in which nexrad files are stored
    """
    
    # Create lists to store wspd, u, and v DataArrays
    refl_all = []
    # Create a list to store the timestamps
    times = []
    
    # Loop over files in the directory
    for filename in os.listdir(drad):
        # read in radar file
        radar = pyart.io.read(f"{drad}{filename}")