#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The Purpose of this Script is to Extract Relevant Lidar Data Based on the 
Available Daily Files and Return this to Main Scritp for Further Analyses
###############################################################################
Correction was Implemented to Account for the Anomalous Data Availability 
Resulting from Overlapping Range Gates

Created on Fri Oct 25 10:28:50 2019
@author: jduncan
@modified: dmurray - modified for CU CubeWind lidar files (WFIP3)

Loading in Doppler Lidar Information:
    Variables Used:
    ###############
    int64 time(time), float32 height(height),
    float32 um(time,height), float32 vm(time,height), float32 wm(time,height),
    float32 Vhm(time,height), float32 Azim(time,height), 
    ###########################################################################
    float32 lat(), float32 lon(), float32 alt()

Corrections were Made to Ensure the Outputted Data is Temporally Sorted
"""

import os, warnings
import netCDF4 as netcdf
import numpy as np
import numpy.matlib
from datetime import datetime,timezone,time
#from pslutils import met2cart
from met_to_uv import wd2uv

def dl_ingest(dir_load,avail_flist,doi):
    #input includes the path to files, available file list, and date of interest
    fil2proc = [fname for fname in avail_flist if doi in fname]
    ###########################################################################
    if len(fil2proc) == 0:
        #No Measurement Files for this Data -- Return Empty Tuples
        DLloc = (); DLxy = ();
        DLdata = ();
        ########################
        METloc = (); METxy = ();
        METdata = ();
        #######################################################################
    else:
        bt=datetime.strptime(doi,'%Y%m%d').replace(tzinfo=timezone.utc)
        midnight_time=datetime.timestamp(bt) # midnight of the day secs since 1970-01-01
        if len(fil2proc) == 1:
            #One Daily File Exists
            DLfil = netcdf.Dataset(os.path.join(dir_load,fil2proc[0]))  
            vars=DLfil.variables
            ###################################################################
            #Timing Information
            time = np.asarray(DLfil.variables['time']) #secs since 1970-01-01
            time_meas=time-midnight_time  # base time offset from midnight
            ###################################################################
            #Measurement Height (AGL)
            hgt = np.asarray(DLfil.variables['height']) #measurement height above ground level
            ###################################################################
            #Velocity Information
            Wvel = np.asarray(DLfil.variables['wm']) #vertical component of wind vector
            WS = np.asarray(DLfil.variables['Vhm']) #wind speed
            WD = np.asarray(DLfil.variables['Azim']) #wind direction
            ###################################################################
            #DL Location Information
            DLlon = np.asarray(DLfil.variables['lon']) #East Longitude
            DLlat = np.asarray(DLfil.variables['lat']) #North Latitude
            DLalt = np.asarray(DLfil.variables['alt']) #Altitutde above MSL
            ###################################################################
                
        else:
            #Multiple Daily Files -- Must Concatenate into Single Arrays
            ###################################################################
            # Sort to Ensure Data Files	Correctly Concatenated #
            ftime = [mfil.split('.')[3] for mfil in fil2proc]
            fil2proc = [x for _,x in sorted(zip(np.asarray(ftime),fil2proc))]
            ###################################################################
            for file_run in fil2proc:
                ###############################################################
                DLfil = netcdf.Dataset(os.path.join(dir_load,file_run))  
                ###############################################################
                #Save Static Variables
                if file_run == fil2proc[0]:
                    hgt = np.asarray(DLfil.variables['height'])
                    DLlon = np.asarray(DLfil.variables['lon'])
                    DLlat = np.asarray(DLfil.variables['lat'])
                    DLalt = np.asarray(DLfil.variables['alt']) 
                    ###########################################################
                    #Timing info
                    time = np.asarray(DLfil.variables['time']) #secs since 1970-01-01
                    time_meas=time-midnight_time  # base time offset from midnight
                    ###########################################################
                    #Velocity Information
                    Wvel = np.asarray(DLfil.variables['wm']) #vertical component of wind vector
                    WS = np.asarray(DLfil.variables['Vhm']) #wind speed
                    WD = np.asarray(DLfil.variables['Azim']) #wind direction
                else:
                    #Load in the Remainder of the Information (Concatenate in Time Dimension)
                    time = np.asarray(DLfil.variables['time']) #secs since 1970-01-01
                    tmeas=time-midnight_time  # base time offset from midnight
                    time_meas = np.concatenate((time_meas,tmeas),axis = 0)
                    ###########################################################
                    Wvel = np.concatenate((Wvel,np.asarray(DLfil.variables['wm'])),axis = 0)
                    WS = np.concatenate((WS,np.asarray(DLfil.variables['Vhm'])),axis = 0)
                    WD = np.concatenate((WD,np.asarray(DLfil.variables['Azim'])),axis = 0)
                    ###########################################################

        #######################################################################        
        #Account for Missing Values (Flagged with -9999 Value)       
        #######################################################################
        WS[WS == -9999.0] = np.nan; WD[WD ==-9999.0] = np.nan
        Wvel[Wvel == -9999.0] = np.nan
        #  Measured U/V are in math coords, so we calculate from WS/WD
        Uvel,Vvel = wd2uv(WS,WD);
        #######################################################################
        WS = np.transpose(WS); WD = np.transpose(WD)
        Uvel = np.transpose(Uvel); Vvel = np.transpose(Vvel); Wvel = np.transpose(Wvel)
        #######################################################################
            
        """ Ensure Arrays are Temporally Sorted Prior to Output """
        if not all(time_meas[i] <= time_meas[i+1] for i in range(len(time_meas)-1)):
            # If Not Sorted -- Sort Data Input #
            SORTind = time_meas.argsort(); # define sorted indices
            ###################################################################
            time_meas = time_meas[SORTind]
            time_bounds = time_bounds[SORTind,:]
            ###################################################################
            WS = WS[:,SORTind]; WD = WD[:,SORTind];
            Uvel = Uvel[:,SORTind]; Vvel = Vvel[:,SORTind]; Wvel = Wvel[:,SORTind];
            ###################################################################

        ###################################################################
        #Return the Data of Interest -- But First Create Tuples
        DLloc = (DLlon,DLlat,DLalt); DLxy = (time_meas,hgt); DLdata = (WS,WD,Uvel,Vvel,Wvel)
        METloc = (); METxy = (); METdata = ();
        ###################################################################
            
    return DLloc, DLxy, DLdata, METloc, METxy, METdata
        
   
