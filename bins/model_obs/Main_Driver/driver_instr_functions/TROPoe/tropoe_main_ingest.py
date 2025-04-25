#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The Purpose of this Script is to Extract Relevant Radar Wind Profiler Data Based
on the Available Daily Files and Return this to Main Scritp for Further Analyses
###############################################################################
Currently Analyzing One-Hour Average Data -- Additional TROPoe Available Data 
Includes 'Hi-Res' Data Information. This Datastream Contains the Same Variables
as Those Denoted Below -- Except at a Higher (i.e. 5-min) Temporal Resolution 
Derived from a Moving (10-min Avg) Window

Created on Fri Oct 25 10:28:50 2019
@author: jduncan

Note that this Ingest Script Currently Assumes that if there are Multiple Daily
Files that they Have the Same Measurement Heights. This Script Only Extracts 
The Measurement Heights from the First Daily File Analyzed

Loading in Radar Wind Profiler Information Information:
    Variables Used:
    ###############
    int base_time, double time_offset(time), 
    float height(height), float temperature(time, height), float waterVapor(time, height)
    float pressure(time, height) float theta(time, height) float thetae(time, height) float rh(time, height) float dewpt(time, height)
    float lwp(time) float gamma(time) float rmsr(time) float pblh(time) float pwv(time) float cbh(time)
    ###########################################################################
    float lat, float lon, float alt
    
    Variables Not Currenlty Used:
    #############################
    float lReff(time) float iTau(time) float iReff(time) float co2(time, gas_dim)
    float ch4(time, gas_dim) float n2o(time, gas_dim) float sigma_temperature(time, height) float sigma_waterVapor(time, height)
    float sigma_lwp(time) float sigma_lReff(time) float sigma_iTau(time) float sigma_iReff(time) float sigma_co2(time, gas_dim)
    float sigma_ch4(time, gas_dim) float sigma_n2o(time, gas_dim)  float rmsa(time) float rmsp(time)
    float chi2(time) float convergence_criteria(time) float dfs(time, dfs) float sic(time) float vres_temperature(time, height)
    float vres_waterVapor(time, height) float cdfs_temperature(time, height) float cdfs_waterVapor(time, height) 
    float dindices(time, index_dim) float sigma_dindices(time, index_dim) float obs_vector(time, obs_dim) float obs_vector_uncertainty(time, obs_dim)
    float forward_calc(time, obs_dim) float Xop(time, arb) float Akernal(time, arb, arb) float Xa(arb)

Corrections:
- Corrections were Made to Ensure the Outputted Data is Temporally Sorted
- Corrections were Made to Only Incorporate Data if Two Scanning Modes Were Employed
"""

import os, warnings
import netCDF4 as netcdf
import numpy as np
import numpy.matlib
from datetime import datetime, timezone, timedelta

def tropoe_ingest(dir_load,avail_flist,doi,suffix='.nc'):
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
    fil_lev_num = [nc.dimensions['height'].size for nc in fil_nc];
    # previously required exactly 121 levels -- but these vary depending on settings
    ###########################################################################
    hgt_chck = len(np.unique(np.asarray(fil_lev_num))) > 1;
    # indicates whether internal concatenation of files can be performed
    ###########################################################################
    # denote reference for number of levels to only concatenate correct data
    if len(fil_lev_num)>0:
        lev2ref = fil_lev_num[0];
    else:
        lev2ref = 0;
    ###########################################################################

    
    if len(fil2proc) == 0 or hgt_chck:
        #No Measurement Files for this Data -- Return Empty Tuples
        TROPoeloc = (); TROPoexy = ()
        TROPoetemp_rh = (); TROPoepres = (); TROPoetheta = (); 
        TROPoemisc = ();
    else:
        #If Datastreams Exist -- Process
        #If multiple files for same day, select the latest version
        if len(fil2proc) > 1:
           fil2proc.sort(reverse=True)
        #Select the first file in the list
        TROPoefil = netcdf.Dataset(os.path.join(dir_load,fil2proc[0]))  
        vars=TROPoefil.variables
        #######################################################################
        #Timing Information
        ###################################################################
        if ('time' in vars) and False:  # older files don't have time variable
            time_meas = np.asarray(TROPoefil.variables['time']) #time offset from midnight (i.e sec since 00:00:00 0:00)
        else:   #create base_time and time_offset from datetime
            bt=datetime.strptime(doi,'%Y%m%d').replace(tzinfo=timezone.utc)
            midnight_time=datetime.timestamp(bt) # midnight of the day secs since 1970-01-01
            base_time = np.asarray(TROPoefil.variables['base_time']) #secs since 1970-01-01
            offset_from_midnight=base_time-midnight_time  # base time offset from midnight
            time_offset = np.asarray(TROPoefil.variables['time_offset']) #time offset from base_time
            time_meas = time_offset+offset_from_midnight

        #Measurement Height (AGL)
        hgt = np.asarray(TROPoefil.variables['height']) #measurement height above ground level (km)
        ###################################################################
        # Denote Wind Measurement Quality #
        Dqual = np.asarray(TROPoefil.variables['qc_flag']) # 0 is okay, threshold of 6
        ###################################################################
        temp = np.asarray(TROPoefil.variables['temperature']) #measurement height above ground level (km)
        rh = np.asarray(TROPoefil.variables['rh']) #measurement height above ground level (km)
        dewpt = np.asarray(TROPoefil.variables['dewpt']) #measurement height above ground level (km)
        wv = np.asarray(TROPoefil.variables['waterVapor']) #measurement height above ground level (km)
        ###################################################################
        theta = np.asarray(TROPoefil.variables['theta']) #measurement height above ground level (km)
        thetae = np.asarray(TROPoefil.variables['thetae']) #measurement height above ground level (km)
        ###################################################################
        pres = np.asarray(TROPoefil.variables['pressure']) #measurement height above ground level (km)
        ###################################################################
        #  Misc vars
        cbh = np.asarray(TROPoefil.variables['cbh']) #cloud base height above ground level (km)
        gamma = np.asarray(TROPoefil.variables['gamma']) # gamma parameter
        rmsr = np.asarray(TROPoefil.variables['rmsr']) # Root mean square error between IRS and MWR obs in the observation vector and the forward calculation
        lwp = np.asarray(TROPoefil.variables['lwp']) # liquid water path
        # Misc Information:  NB: for older files, these are stored in dindices
        #precipitable water
        if ('pwv' in vars):
          pwv = np.asarray(TROPoefil.variables['pwv']); 
        else:  
          #see if it's in the dindices var
          if ('dindices' in vars):
            temp_ = np.asarray(TROPoefil.variables['dindices'])
            pwv = temp_[:,0] #precipitable water
            del temp_
          else:
            pwv = np.full(temp.shape,np.nan);
        #planetary boundary layer
        if ('pblh' in vars):
          pblh = np.asarray(TROPoefil.variables['pblh']); 
        else:
          #see if it's in the dindices var
          if ('dindices' in vars):
            temp_ = np.asarray(TROPoefil.variables['dindices'])
            pblh = temp_[:,1]
            del temp_
          else:
            pblh = np.full(temp.shape,np.nan);

        ###################################################################
        TROPoelon = np.asarray(TROPoefil.variables['lon']) #East Longitude
        TROPoelat = np.asarray(TROPoefil.variables['lat']) #North Latitude
        TROPoealt = np.asarray(TROPoefil.variables['alt']) #Altitutde above MSL            
        
        #######################################################################
        #Account for Missing Values (Flagged with -9999 Value)
        #######################################################################
        hgt[hgt == 999999.0] = np.nan;
        hgt = hgt*1000 # convert to meters
        #######################################################################
        temp[temp == 999999.0] = np.nan; rh[rh == 999999.0] = np.nan;
        dewpt[dewpt == 999999.0] = np.nan; wv[wv == 999999.0] = np.nan;
        theta[theta == 999999.0] = np.nan; thetae[thetae == 999999.0] = np.nan;
        pres[pres == 999999.0] = np.nan; 
        #######################################################################
        # Account for the way python interprets matrices (i.e. row x column)
        temp = np.transpose(temp); rh = np.transpose(rh)
        dewpt = np.transpose(dewpt); wv = np.transpose(wv)
        theta = np.transpose(theta); thetae = np.transpose(rh)
        pres = np.transpose(pres); 
    
        """ Package Tuples for Output """
        TROPoeloc = (TROPoelon,TROPoelat,TROPoealt); TROPoexy = (time_meas,hgt);
        #######################################################################
        #TROPoedata = (WS,WD,Uvel,Vvel,Wvel,SNRbeam,SWbeam,Wqual)
        TROPoetemp_rh = (temp,rh,dewpt,wv,pres)
        TROPoetheta = (theta,thetae)
        TROPoemisc = (cbh,gamma,rmsr,lwp,pwv,pblh)

    return TROPoeloc, TROPoexy, TROPoetemp_rh, TROPoetheta, TROPoemisc, lev2ref
