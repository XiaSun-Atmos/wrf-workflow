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
    >f8 time_offset(time)    (seconds since midnight)
    >f4 pressure(time)
    >f4 temperature(time) 
    >f4 relative_humidity(time)
    >f4 wind_speed(time)
    >f4 wind_direction(time)
    >f4 solar_radiation(time)
    >f4 net_radiation(time)
    >f4 lat(), >f4 lon(), >f4 alt()
    
"""

import os
import netCDF4 as netcdf
import numpy as np

def pslmet_ingest(dir_load,avail_flist,doi):
    #input includes the path to files, available file list, and date of interest
    fil2proc = [fname for fname in avail_flist if ((doi in fname) & (fname.endswith('.nc')))]
    ###########################################################################
    #fil_nc = [netcdf.Dataset(os.path.join(dir_load,f)) for f in fil2proc]
    fil_nc = []
    for f in fil2proc:
       try:
          fil_nc.append(netcdf.Dataset(os.path.join(dir_load,f)))
       except:
          print('unable to open',f)
          fil2proc.remove(f)
          continue
    ###########################################################################
    if len(fil2proc) == 0:
        #No Measurement Files for this Data -- Return Empty Tuples
        PSLloc = (); PSLxy = ()
        PSLwind = (); PSLtemp_rh = (); PSLpres = ()
        #######################################################################
        PSLprecip = ()
        #######################################################################
        PSLrad = ()
    else:
        #If Datastreams Exist -- Process
        if len(fil2proc) == 1:
            #One Daily File Exists
            PSLfil = netcdf.Dataset(os.path.join(dir_load,fil2proc[0]))  
            vars=PSLfil.variables
            #######################################################################
            #Timing Information
            time_meas = np.asarray(PSLfil.variables['time_offset']) #time offset from midnight (i.e sec since 00:00:00 0:00)
            
            #Atmo Information Below:
            #######################################################################
            #Pressure Information:
            ATMOpress = np.asarray(PSLfil.variables['pressure']);
            #######################################################################
            #Temperature Information:
            ATMOtemp = np.asarray(PSLfil.variables['temperature']);
            #######################################################################
            #Relative Humidity Information:
            ATMOrh = np.asarray(PSLfil.variables['relative_humidity']);
            #######################################################################
            #Wind Speed Information:
            if ('wind_speed' in vars):
              ATMOws = np.asarray(PSLfil.variables['wind_speed']);
            else:
              ATMOws = np.full(ATMOtemp.shape,np.nan);
            #######################################################################
            #Wind Direction Information:
            if ('wind_direction' in vars):
              ATMOwd_vec = np.asarray(PSLfil.variables['wind_direction']);
            else:
              ATMOwd_vec = np.full(ATMOtemp.shape,np.nan);
            if ('stdwind_direction' in vars):
              ATMOwd_vec_std = np.asarray(PSLfil.variables['stdwind_direction']);
            else:
              ATMOwd_vec_std = np.full(ATMOtemp.shape,np.nan);
            #######################################################################
            #Tipping Bucket Raing Gauge Precipitation Information:
            if ('precipitation' in vars):
              ATMOprecipTBRG_total = np.asarray(PSLfil.variables['precipitation'])
            else:
              ATMOprecipTBRG_total = np.full(ATMOtemp.shape,np.nan);
            #######################################################################
            #Radiation Information:
            if ('solar_radiation' in vars):
              ATMOsrad = np.asarray(PSLfil.variables['solar_radiation']);
            else:
              ATMOsrad = np.full(ATMOtemp.shape,np.nan);
            if ('net_radiation' in vars):
              ATMOnrad = np.asarray(PSLfil.variables['net_radiation']);
            else:
              ATMOnrad = np.full(ATMOtemp.shape,np.nan);
            #######################################################################
            #PSL Location Information
            PSLlon = np.asarray(PSLfil.variables['lon']) #East Longitude
            PSLlat = np.asarray(PSLfil.variables['lat']) #North Latitude
            PSLalt = np.asarray(PSLfil.variables['alt']) #Altitutde above MSL
  
        #######################################################################    
        #Account for Missing Values (Flagged with -9999 Value)    
        #-> Do Not Apply to QC Variables 
        #######################################################################
        ATMOpress[ATMOpress == 99999.0] = np.nan; 
        ATMOtemp[ATMOtemp == 99999.0] = np.nan; 
        ATMOrh[ATMOrh == 99999.0] = np.nan; 
        #######################################################################
        ATMOws[ATMOws == 99999.0] = np.nan; 
        ATMOwd_vec[ATMOwd_vec == 99999.0] = np.nan; 
        ATMOwd_vec_std[ATMOwd_vec_std == 99999.0] = np.nan;
        #######################################################################
        #Data Analyses Indicates that when the Arithmetic Wind Speed is Equal to 0,
        #Then the Wind Direction is Also Set to Zero (Do Not Consider These Obs)
        ATMOwd_vec[ATMOws == 0] = np.nan; ATMOwd_vec_std[ATMOws == 0] = np.nan;
        #######################################################################
        ATMOprecipTBRG_total[ATMOprecipTBRG_total == 99999.0] = np.nan; 
        #######################################################################
        ATMOsrad[ATMOsrad == 99999.0] = np.nan; 
        ATMOnrad[ATMOnrad == 99999.0] = np.nan; 
        #######################################################################
        #Return the Data of Interest -- But First Create Tuples
        PSLloc = (PSLlon,PSLlat,PSLalt)
        #-> Create Tuple by Atmosphereic Characterization
        PSLxy = (time_meas)
        PSLwind = (ATMOws,ATMOwd_vec,ATMOwd_vec_std)
        PSLtemp_rh = (ATMOtemp,ATMOrh)
        PSLpres = (ATMOpress)
        #######################################################################
        TBRGprecip = (ATMOprecipTBRG_total)
        PSLprecip = (TBRGprecip) 
        PSLrad = (ATMOsrad,ATMOnrad)
        #######################################################################
    
    return PSLloc, PSLxy, PSLwind, PSLtemp_rh, PSLpres, PSLprecip, PSLrad
