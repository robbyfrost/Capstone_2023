# ------------------------------------------------
# Name: run_utils.py
# Author: Robert M. Frost, Savannah J. Southward, and Colin K. Welty
# University of Oklahoma
# Created: 09 November 2023
# Purpose: Script to run functions from 
# capstone_utils.py
# ------------------------------------------------
from capstone_utils import *

# settings
run_cartesian_column = False
run_timeheight = True

drad = "/Users/robertfrost/Documents/boomsoon/2023_semester_2/capstone/analysis/Capstone_2023/data/lightning_cases/201708-055/"
dtrack = "/Users/robertfrost/Documents/boomsoon/2023_semester_2/capstone/analysis/capstone_2023.nc"
meso_id = "201708-055"
dout = "/Users/robertfrost/Documents/boomsoon/2023_semester_2/capstone/analysis/Capstone_2023/data/lightning_cases/"
track_point = 0

# run functions

# extract column for single time
if run_cartesian_column:
    cartesian_column(drad, dtrack, meso_id, track_point)

# create time-height series
if run_timeheight:
    timeheight(drad, dtrack, meso_id, dout)