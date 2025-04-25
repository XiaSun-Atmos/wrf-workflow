#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The Purpose of this Function is to Properly Concatenate Input Variables and to
Return the Array of Interest

Inputted into the function is the Tuple Information for the Day Before, the 
Current Day, and the Following Day
-- An empty tuple is inherently inputted if that data file did not exist

** Input Tuples Must Correclty Correspond to the Day Being Analyzed (i.e. Pre,
** Day, and Pst)

Script was Updated to Handle 'Extended Forecast Periods' That Require an EXT 
Data File to be Processed
###############################################################################
Created on Tue Dec 31 14:25:52 2019
@author: jduncan
"""
##################
import numpy as np

def atmo_cat(pre,day,pst,ext,var_oi):
    """ Properly Concatenate Input Based On Available Information """
    INPtuple_arr = [pre,day,pst,ext] # input time arrays
    ###########################################################################
    CATarr = [];
    for ioi in range(0,len(INPtuple_arr)):
        if len(INPtuple_arr[ioi]) > 0:
            if var_oi >= 0:
              CATarr.append(INPtuple_arr[ioi][var_oi])
            else:
              CATarr.append(INPtuple_arr[ioi])
    #CATarr = [x for x in INPtuple if len(x) > 0]
    ###########################################################################
    if len(CATarr) == 0:
        ATMOanly = ()
    elif len(CATarr) == 1:
        ATMOanly = CATarr[0]
    elif len(CATarr) == 2:
        ATMOanly = np.concatenate((CATarr[0],CATarr[1]),axis=0)
    elif len(CATarr) == 3:
        ATMOanly = np.concatenate((CATarr[0],CATarr[1],CATarr[2]),axis=0)
    else:
        ATMOanly = np.concatenate((CATarr[0],CATarr[1],CATarr[2],CATarr[3]),axis=0)
    # the efficiency of this code could be improved
    ###########################################################################
    return ATMOanly
