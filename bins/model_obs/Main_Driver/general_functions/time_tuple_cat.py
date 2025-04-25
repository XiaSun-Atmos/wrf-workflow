#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The Purpose of this Function is to Properly Concatenate Time Input Variables 
and to Return the Concatenated Time Array

Inputted into the function is the Tuple Information for the Day Before, the 
Current Day, and the Following Day
-- An empty tuple is inherently inputted if that data file did not exist
###############################################################################
This Function Can Handle:
    - Tuple Input -- i.e. the time is also combined with time_bounds or other 
      similar variables. Require length of tuple to be greater than 0 to 
      proceed given a 0 input tuple indicates no data
    - Represents Cases Where the Only Data Provided in the Initial xy Tuple was 
      the time_meas data

Script was Updated to Handle 'Extended Forecast Periods' That Require an EXT 
Data File to be Processed
###############################################################################
Created on Tue Dec 31 14:25:52 2019
@author: jduncan
"""
##################
import numpy as np
##################
HRsec = 60*60;
DAYsec = 24*HRsec;
EXTsec = DAYsec*2;

def time_cat(pre,day,pst,ext):
    """ Properly Concatenate Input Based On Available Information """
    INPtuple_arr = [pre,day,pst,ext] # input time arrays
    INPtuple_off = [-DAYsec,0,DAYsec,EXTsec] # required time offset
    ###########################################################################
    CAToff = [];
    CATarr = [];
    ###########################################################################
    for ioi in range(0,len(INPtuple_arr)):
        #######################################################################
        # Determine Whether Item is Numpy Array (Meaning Only Time Data was 
        # Provided in the pre,day,pst arrays)
        if isinstance(INPtuple_arr[ioi],np.ndarray):
            CATarr.append(INPtuple_arr[ioi])
            CAToff.append(INPtuple_off[ioi])    
        elif isinstance(INPtuple_arr[0],tuple) and len(INPtuple_arr[ioi]) > 1:
            #indicating it is not an empty tuple -- subsequently reference 0 
            #indice for time data
            CATarr.append(INPtuple_arr[ioi][0])
            CAToff.append(INPtuple_off[ioi])     
    # neglect tuple input with length of 0 (i.e. empty tuple)
    ###########################################################################
    #CATarr = [x for x in INPtuple if len(x) > 0]
    ###########################################################################
    if len(CATarr) == 0:
        TIMEanly = ()
    elif len(CATarr) == 1:
        TIMEanly = CATarr[0]+CAToff[0]
    elif len(CATarr) == 2:
        TIMEanly = np.concatenate((CATarr[0]+CAToff[0],CATarr[1]+CAToff[1]),axis=0)
    elif len(CATarr) == 3:
        TIMEanly = np.concatenate((CATarr[0]+CAToff[0],CATarr[1]+CAToff[1],CATarr[2]+CAToff[2]),axis=0)
    else:
        TIMEanly = np.concatenate((CATarr[0]+CAToff[0],CATarr[1]+CAToff[1],CATarr[2]+CAToff[2],\
                                   CATarr[3]+CAToff[3]),axis=0)
    # the efficiency of this code could be improved
    ###########################################################################
    return TIMEanly
