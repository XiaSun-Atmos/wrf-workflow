#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The Purpose of this Function is to Properly Concatenate Input Variables and to
Return the Array of Interest

Inputted into the function is the Tuple Information for the Day Before, the 
Current Day, and the Following Day
-- An empty tuple is inherently inputted if that data file did not exist

This script differs from atmo_cat in that we are concatenating a 2-d field as 
opposed to a 1-d field
--  To Solve this Issue Simply Concatenate the Data Along the Time Axis Defined
    by the Input

###############################################################################
Created on Tue Dec 31 14:25:52 2019
@author: jduncan
"""
##################
import numpy as np

def atmo_spatial_cat(pre,day,pst,ext,var_oi,time_ind):
    """ Properly Concatenate Input Based On Available Information """
    INPtuple_arr = [pre,day,pst,ext] # input time arrays
    ###########################################################################
    CATarr = [];
    for ioi in range(0,len(INPtuple_arr)):
        if len(INPtuple_arr[ioi]) > 0:
            CATarr.append(INPtuple_arr[ioi][var_oi])
    #CATarr = [x for x in INPtuple if len(x) > 0]
    ###########################################################################
    if len(CATarr) == 0:
        ATMOanly = ()
    elif len(CATarr) == 1:
        ATMOanly = CATarr[0]
    elif len(CATarr) == 2:
        ATMOanly = np.concatenate((CATarr[0],CATarr[1]),axis=time_ind)
    elif len(CATarr) == 3:
        ATMOanly = np.concatenate((CATarr[0],CATarr[1],CATarr[2]),axis=time_ind)
    else:
        ATMOanly = np.concatenate((CATarr[0],CATarr[1],CATarr[2],CATarr[3]),axis=time_ind)
    ###########################################################################
    return ATMOanly
