#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The Purpose of this Script is to Extract Relevant PSL Meteorological Information
Given the Inputted Data File
-> This Script is Able to Handle Multiple Daily Files if Necessary
###############################################################################
Created on Fri Oct 25 10:28:50 2019
@author: jduncan
@modified: dmurray

Loading in PSL Standard Meteorological Information
    Variables Used:
    ###############
    >f8 time(time)
    >f4 Rate(time)
    >f4 Amount(time)
    >f4 AmountSum(time)
    >f4 lat(), >f4 lon(), >f4 alt()
    
"""

import os
import netCDF4 as netcdf
import numpy as np

def psldisd_ingest(dir_load,avail_flist,doi):
    #input includes the path to files, available file list, and date of interest
    fil2proc = [fname for fname in avail_flist if doi in fname]
    ###########################################################################
    if len(fil2proc) == 0:
        #No Measurement Files for this Data -- Return Empty Tuples
        PSLloc = (); PSLxy = ()
        #######################################################################
        PSLprecip = ()
    else:
        PSLprecip = ()
        #If Datastreams Exist -- Process
        if len(fil2proc) == 1:
            #One Daily File Exists
            PSLfil = netcdf.Dataset(os.path.join(dir_load,fil2proc[0]))  
            #######################################################################
            #Timing Information
            time_meas = np.asarray(PSLfil.variables['time_offset']) #time offset from midnight (i.e sec since 00:00:00 0:00)
            time_meas_start = np.asarray(PSLfil.variables['time']) #start time of measurement period
            time_meas_end = np.asarray(PSLfil.variables['time_end']) #end time of measurement period
            time_meas_interval = np.asarray(PSLfil.variables['time_interval']) #end time of measurement period
            
            #Atmo Information Below:
            #######################################################################
            #Precip Rate Information:
            ATMOprate = np.asarray(PSLfil.variables['Rate']);
            #long_name: Precipitation Rate
            #units: mm/hr
            #######################################################################
            #Precip Amount Information:
            ATMOpr = np.asarray(PSLfil.variables['Amount']);
            #long_name: Precipitation Amount
            #units: mm
            #######################################################################
            #Relative Humidity Information:
            ATMOprsum = np.asarray(PSLfil.variables['AmountSum']);
            #long_name: Event Precipitation Amount 
            #units: mm
            #######################################################################
            #PSL Location Information
            PSLlon = np.asarray(PSLfil.variables['lon']) #East Longitude
            PSLlat = np.asarray(PSLfil.variables['lat']) #North Latitude
            PSLalt = np.asarray(PSLfil.variables['alt']) #Altitutde above MSL
  
        #######################################################################    
        #Account for Missing Values (Flagged with -9999 Value)    
        #-> Do Not Apply to QC Variables 
        #######################################################################
        ATMOprate[ATMOprate == 99999] = np.nan; 
        ATMOpr[ATMOpr == 99999] = np.nan; 
        ATMOprsum[ATMOprsum == 99999] = np.nan; 
        #######################################################################
        #Return the Data of Interest -- But First Create Tuples
        PSLloc = (PSLlon,PSLlat,PSLalt)
        #-> Create Tuple by Atmosphereic Characterization
        PSLxy = (time_meas,time_meas_start,time_meas_end,time_meas_interval)
        PSLprecip = (ATMOprate,ATMOpr,ATMOprsum)
        #######################################################################
        #######################################################################
    
    return PSLloc, PSLxy, PSLprecip
