#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The Purpose of this Script is to Extract GML Ceilometer data
###############################################################################
Created on Tue, Feb 20, 2024
@author: dmurray

Loading in GML Ceilometer Information:
    Variables Used:
    ###############
    int base_time, double time_offset(time), 
    float first_cbh(time), 
"""
    
import os, warnings
from datetime import datetime
import netCDF4 as netcdf
import numpy as np
import numpy.matlib

def gmlceil_ingest(dir_load,avail_flist,doi,suffix='.cdf'):
    #input includes the path to files, available file list, and date of interest
    fil2proc = [fname for fname in avail_flist if ((doi in fname) & (fname.endswith(suffix)))]
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
    
    if len(fil2proc) == 0:
        #No Measurement Files for this Data -- Return Empty Tuples
        CEILloc = (); CEILxy = ()
        CEILhgt = (); 
    else:
        #If Datastreams Exist -- Process
        #If multiple files for same day, select the latest version
        if len(fil2proc) > 1:
           fil2proc.sort(reverse=True)
        #Select the first file in the list
        CEILfil = netcdf.Dataset(os.path.join(dir_load,fil2proc[0]))  
        vars=CEILfil.variables
        #######################################################################
        #Timing Information
        ###################################################################
        time_base = np.asarray(CEILfil.variables['base_time']) #time offset from midnight (i.e sec since 00:00:00 0:00)
        time_meas = np.asarray(CEILfil.variables['time']) #time offset from midnight (i.e sec since 00:00:00 0:00)
        #time_meas=np.asarray([time_base+time_offset[i] for i in range(len(time_offset))])

        #Measurement Height (AGL)
        hgt = np.asarray(CEILfil.variables['first_cbh']).astype(np.float32) #measurement height above ground level (m)
        hgt[hgt < 0] = np.nan;

        """ Package Tuples for Output """
        CEILloc = (); CEILxy = (time_meas,hgt);
        #######################################################################
        CEILhgt = (hgt)

    return CEILloc, CEILxy, CEILhgt
