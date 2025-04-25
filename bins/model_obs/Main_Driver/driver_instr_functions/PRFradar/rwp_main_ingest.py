#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The Purpose of this Script is to Extract Relevant Radar Wind Profiler Data Based
on the Available Daily Files and Return this to Main Scritp for Further Analyses
###############################################################################
Currently Analyzing One-Hour Average Data -- Additional RWP Available Data 
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
    float64 time(time), float64 time_bounds(time,bound). float32 height(mode,level), 
    float32 wind_speed(time,mode,level), float32 wind_speed_std(time,mode,level), 
    float32 wind_direction(time,mode,level), float32 wind_quality(time,mode,level), 
    float32 wind_u(time,mode,level), float32 wind_v(time,mode,level), float32 wind_w(time,mode,level), 
    float32 radial_velocity_snr(time,mode,beam,level),float32 radial_velocity_width(time,mode,beam,level), 
    ###########################################################################
    float32 lat(), float32 lon(), float32 alt()   
    
    Variables Not Currenlty Used:
    #############################
    int32 base_time(), float64 time_offset(time), S1 mode_name(mode,strlen), int32 mode_number(mode),   
    int32 transmit_power_level_setting(mode), float32 pulse_width(mode), float32 interpulse_period(mode),
    int32 fourier_transform_length(mode), int32 num_time_domain(mode), int32 num_frequency_domain(mode),
    int32 averaging_interval(mode), int32 wind_quality_interval(mode), int32 num_range_gates(mode), 
    int32 num_levels(mode), int32 num_beams(mode), float32 wind_w_std(time,mode,level), 
    float32 zenith_angle(mode), float32 azimuth(mode,beam),float32 radial_velocity(time,mode,beam,level), 
    int32 radial_velocity_navg(time,mode,beam,level),float32 radial_velocity_power(time,mode,beam,level), 

Corrections:
- Corrections were Made to Ensure the Outputted Data is Temporally Sorted
- Corrections were Made to Only Incorporate Data if Two Scanning Modes Were Employed
"""

import os, warnings
import netCDF4 as netcdf
import numpy as np
import numpy.matlib

def rwp_ingest(dir_load,avail_flist,doi):
    #input includes the path to files, available file list, and date of interest
    fil2proc = [fname for fname in avail_flist if doi in fname]
    ###########################################################################
    # Adjust fil2proc to Ensure Each Explictly Has Two Scanning Modes
    fil_nc = [netcdf.Dataset(os.path.join(dir_load,f)) for f in fil2proc]
    fil_mode_num = [nc.dimensions['mode'].size for nc in fil_nc];
    fil_lev_num = [nc.dimensions['level'].size for nc in fil_nc];
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
    fil2proc = [fil2proc[f] for f in range(len(fil2proc)) if fil_mode_num[f] == 2]
    # using this method should preserve those files with two scanning modes #
    
    if len(fil2proc) == 0 or hgt_chck:
        #No Measurement Files for this Data -- Return Empty Tuples
        RWPloc = (); RWPxy = ();
        RWPdata = (); RWPlow = (); RWPhgh = (); lev2ref = 0;

    else:
        if len(fil2proc) == 1:
            #One Daily File Exists
            RWPfil = netcdf.Dataset(os.path.join(dir_load,fil2proc[0]))  
            ###################################################################
            # Define the Number of Range Gates Associated with the Different
            # Sampling Modes
            LOWmode = np.asarray(RWPfil.variables['num_range_gates'][0])
            HGHmode = np.asarray(RWPfil.variables['num_range_gates'][1])
            ###################################################################
            #Timing Information
            #time_meas = np.asarray(RWPfil.variables['time']) #time offset from midnight (i.e sec since 00:00:00 0:00)
            time_meas = np.asarray(RWPfil.variables['time_offset']) #time offset from midnight (i.e sec since 00:00:00 0:00)
            #time_bounds = np.asarray(RWPfil.variables['time_bounds']) #bounds of the scan times (i.e. start/finish) [i.e. - 3600 s/0 s]
            #Measurement Height (AGL)
            hgt = np.asarray(RWPfil.variables['height']) #measurement height above ground level (km)
            ###################################################################
            # Denote Wind Measurement Quality #
            Wqual = np.asarray(RWPfil.variables['wind_quality']) #Decimal value between 0 and 1 indicating the quality of the winds at the given level
            ###################################################################
            #Velocity Information
            Uvel = np.asarray(RWPfil.variables['wind_u']) #zonal wind speed
            Vvel = np.asarray(RWPfil.variables['wind_v']) #meridonal wind speed
            #Wvel = np.asarray(RWPfil.variables['wind_w']) #vertical wind speed
            ###################################################################
            WS = np.asarray(RWPfil.variables['wind_speed']) #wind speed (time,mode,level) [24,2,121]
            WD = np.asarray(RWPfil.variables['wind_direction']) #wind direction (time,mode,level) [24,2,121]
            ###################################################################
            #Signal-to-Noise Ratio Information
            SNRbeam = np.asarray(RWPfil.variables['radial_velocity_snr']) #Signal to noise ratio of the atmospheric signal associated with radial_velocity
            #Spectral Width
            #SWbeam = np.asarray(RWPfil.variables['radial_velocity_width']) #Spectral width of the atmospheric signal associated with radial_velocity
            ###################################################################
            #RWP Location Information
            RWPlon = np.asarray(RWPfil.variables['lon']) #East Longitude
            RWPlat = np.asarray(RWPfil.variables['lat']) #North Latitude
            RWPalt = np.asarray(RWPfil.variables['alt']) #Altitutde above MSL            
            
####        else:
####            #Multiple Daily Files -- Must Concatenate into Single Arrays
####            ###################################################################
####            # Sort to Ensure Data Files	Correctly Concatenated #
####            ftime = [mfil.split('.')[3] for mfil in fil2proc]
####            fil2proc = [x for _,x in sorted(zip(np.asarray(ftime),fil2proc))]
####            ###################################################################
####            for file_run in fil2proc:
####                ###############################################################
####                RWPfil = netcdf.Dataset(os.path.join(dir_load,file_run))
####                ###############################################################
####                # Save Static Variables
####                if file_run == fil2proc[0]:
####                    # Define the Measurement Location
####                    ###########################################################
####                    RWPlon = np.asarray(RWPfil.variables['lon']) #East Longitude
####                    RWPlat = np.asarray(RWPfil.variables['lat']) #North Latitude
####                    RWPalt = np.asarray(RWPfil.variables['alt']) #Altitutde above MSL
####                    ###########################################################
####                    # Define the Number of Range Gates Assocaited with the Diff
####                    # Sampling Modes
####                    LOWmode = np.asarray(RWPfil.variables['num_range_gates'][0])
####                    HGHmode = np.asarray(RWPfil.variables['num_range_gates'][1])
####                    ###########################################################
####                    time_meas = np.asarray(RWPfil.variables['time']) #time offset from midnight (i.e sec since 00:00:00 0:00)
####                    time_bounds = np.asarray(RWPfil.variables['time_bounds']) #bounds of the scan times (i.e. start/finish) [i.e. - 3600 s/0 s]
####                    hgt = np.asarray(RWPfil.variables['height']) #measurement height above ground level (m)
####                    ###########################################################
####                    Wqual = np.asarray(RWPfil.variables['wind_quality']) #Decimal value between 0 and 1 indicating the quality of the winds at the given level
####                    ###########################################################
####                    Uvel = np.asarray(RWPfil.variables['wind_u']) #zonal wind speed
####                    Vvel = np.asarray(RWPfil.variables['wind_v']) #meridonal wind speed
####                    Wvel = np.asarray(RWPfil.variables['wind_w']) #vertical wind speed
####                    ###########################################################
####                    WS = np.asarray(RWPfil.variables['wind_speed']) #wind speed (time,mode,level) [24,2,121]
####                    WD = np.asarray(RWPfil.variables['wind_direction']) #wind direction (time,mode,level) [24,2,121]
####                    ###########################################################
####                    SNRbeam = np.asarray(RWPfil.variables['radial_velocity_snr']) #Signal to noise ratio of the atmospheric signal associated with radial_velocity
####                    SWbeam = np.asarray(RWPfil.variables['radial_velocity_width']) #Spectral width of the atmospheric signal associated with radial_velocity
####                    ###########################################################
####                else:
####                    #Load in the Remainder of the Information (Concatenate in Time Dimension)
####                    time_meas = np.concatenate((time_meas,np.asarray(RWPfil.variables['time'])),axis = 0)
####                    time_bounds = np.concatenate((time_bounds,np.asarray(RWPfil.variables['time_bounds'])),axis = 0)
####                    ###########################################################
####                    Wqual = np.concatenate((Wqual,np.asarray(RWPfil.variables['wind_quality'])),axis = 0);
####                    ###########################################################
####                    Uvel = np.concatenate((Uvel,np.asarray(RWPfil.variables['wind_u'])),axis = 0)
####                    Vvel = np.concatenate((Vvel,np.asarray(RWPfil.variables['wind_v'])),axis = 0)
####                    Wvel = np.concatenate((Wvel,np.asarray(RWPfil.variables['wind_w'])),axis = 0)
####                    ###########################################################
####                    WS = np.concatenate((WS,np.asarray(RWPfil.variables['wind_speed'])),axis = 0)
####                    WD = np.concatenate((WD,np.asarray(RWPfil.variables['wind_direction'])),axis = 0)
####                    ###########################################################
####                    SNRbeam = np.concatenate((SNRbeam,np.asarray(RWPfil.variables['radial_velocity_snr'])),axis = 0)
####                    SWbeam = np.concatenate((SWbeam,np.asarray(RWPfil.variables['radial_velocity_width'])),axis = 0)
        
        #######################################################################
        #Account for Missing Values (Flagged with -9999 Value)
        #######################################################################
        hgt[hgt == 999999.0] = np.nan;
        hgt = hgt*1000 # convert to meters
        #######################################################################
        #Wqual[Wqual == 999999.0] = np.nan;
        #######################################################################
        Uvel[Uvel == 999999.0] = np.nan; Vvel[Vvel == 999999.0] = np.nan;
        #Wvel[Wvel == 999999.0] = np.nan;
        #######################################################################
        WS[WS == 999999.0] = np.nan; WD[WD == 999999.0] = np.nan;
        #SNRbeam[SNRbeam == 999999.0] = np.nan; SWbeam[SWbeam == 999999.0] = np.nan;
        SNRbeam[SNRbeam == 999999.0] = np.nan; 
        #######################################################################
        # Account for the way python interprets matrices (i.e. row x column)

        """ Prior to Transposing Data -- to Accomodate Python Interpretation --
        Parse Wind Information into the Correspondibng Measurement (i.e. high/
        low altitude) Modes (rows x columns) """
        #######################################################################
        # High- and Low-Altitude Modes #
        LOWind = 0; HGHind = 1;
        # Also Factor in the Relevant Number of Range Gates
        #######################################################################
        HGTlow = hgt[LOWind,0:LOWmode]; HGThgh = hgt[HGHind,0:HGHmode]
        #######################################################################
        Wlow_qual = Wqual[:,LOWind,0:LOWmode];
        Whgh_qual = Wqual[:,HGHind,0:HGHmode];
        ##############################################
        WSlow = WS[:,LOWind,0:LOWmode]; WShgh = WS[:,HGHind,0:HGHmode];
        WDlow = WD[:,LOWind,0:LOWmode]; WDhgh = WD[:,HGHind,0:HGHmode];
        #######################################################################
        Uvel_low = Uvel[:,LOWind,0:LOWmode]; Uvel_hgh = Uvel[:,HGHind,0:HGHmode];
        Vvel_low = Vvel[:,LOWind,0:LOWmode]; Vvel_hgh = Vvel[:,HGHind,0:HGHmode];
        #Wvel_low = Wvel[:,LOWind,0:LOWmode]; Wvel_hgh = Wvel[:,HGHind,0:HGHmode];
        #######################################################################
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            ###################################################################
            SNRlow = np.nanmean(SNRbeam[:,LOWind,:,0:LOWmode],axis=1);
            SNRhgh = np.nanmean(SNRbeam[:,HGHind,:,0:HGHmode],axis=1);
            ###################################################################
            # Peform Range-Correction on SNR Values #
            hgt_lowin = np.matlib.repmat(HGTlow,SNRlow.shape[0],1)
            hgt_hghin = np.matlib.repmat(HGThgh,SNRhgh.shape[0],1)
            # previously height was converted to km for ingest, which is wrong
            # i.e. np.matlib.repmat(HGTlow/1000,SNRlow.shape[0],1)
            # i.e. np.matlib.repmat(HGThgh/1000,SNRhgh.shape[0],1)
            ###################################################################
            SNRlow_rngecorr = SNRlow + 20 * np.log10(hgt_lowin)
            SNRhgh_rngecorr = SNRhgh + 20 * np.log10(hgt_hghin)
            ###################################################################
            #SWlow = np.nanmean(SWbeam[:,LOWind,:,0:LOWmode],axis=1);
            #SWhgh = np.nanmean(SWbeam[:,HGHind,:,0:HGHmode],axis=1);

        #######################################################################
        # Transpose Variables to Conform with Python Interpretation
        Wlow_qual = np.transpose(Wlow_qual); Whgh_qual = np.transpose(Whgh_qual);

        WSlow = np.transpose(WSlow); WShgh = np.transpose(WShgh);
        WDlow = np.transpose(WDlow); WDhgh = np.transpose(WDhgh);

        Uvel_low = np.transpose(Uvel_low); Uvel_hgh = np.transpose(Uvel_hgh);
        Vvel_low = np.transpose(Vvel_low); Vvel_hgh = np.transpose(Vvel_hgh);
        #Wvel_low = np.transpose(Wvel_low); Wvel_hgh = np.transpose(Wvel_hgh);

        SNRlow = np.transpose(SNRlow); SNRhgh = np.transpose(SNRhgh);
        SNRlow_rngecorr = np.transpose(SNRlow_rngecorr); SNRhgh_rngecorr = np.transpose(SNRhgh_rngecorr);
        #SWlow = np.transpose(SWlow); SWhgh = np.transpose(SWhgh);
        
        """ Ensure Arrays are Temporally Sorted Prior to Output """
        if not all(time_meas[i] <= time_meas[i+1] for i in range(len(time_meas)-1)):
            # If Not Sorted -- Sort Data Input #
            SORTind = time_meas.argsort(); # define sorted indices
            ###################################################################
            time_meas = time_meas[SORTind]
            #time_bounds = time_bounds[SORTind,:]
            ###################################################################
            WS = WS[SORTind,:,:]; WD = WD[SORTind,:,:];
            Uvel = Uvel[SORTind,:,:]; Vvel = Vvel[SORTind,:,:]; #Wvel = Wvel[SORTind,:,:];
            Wqual = Wqual[SORTind,:,:];
            ###################################################################
            SNRbeam = SNRbeam[SORTind,:,:,:]; #SWbeam = SWbeam[SORTind,:,:,:];
            ###################################################################
            WSlow = WSlow[:,SORTind]; WDlow = WDlow[:,SORTind];
            Uvel_low = Uvel_low[:,SORTind]; Vvel_low = Vvel_low[:,SORTind]; #Wvel_low = Wvel_low[:,SORTind];
            Wlow_qual = Wlow_qual[:,SORTind];
            ###################################################################
            SNRlow = SNRlow[:,SORTind]; SNRlow_rngecorr = SNRlow_rngecorr[:,SORTind];
            #SWlow = SWlow[:,SORTind];
            ###################################################################
            WShgh = WShgh[:,SORTind]; WDhgh = WDhgh[:,SORTind];
            Uvel_hgh = Uvel_hgh[:,SORTind]; Vvel_hgh = Vvel_hgh[:,SORTind]; #Wvel_hgh = Wvel_hgh[:,SORTind];
            Whgh_qual = Whgh_qual[:,SORTind];
            ###################################################################
            SNRhgh = SNRhgh[:,SORTind]; SNRhgh_rngecorr = SNRhgh_rngecorr[:,SORTind];
            #SWhgh = SWhgh[:,SORTind];

        """ Package Tuples for Output """
        #RWPloc = (RWPlon,RWPlat,RWPalt); RWPxy = (time_meas,time_bounds,hgt);
        RWPloc = (RWPlon,RWPlat,RWPalt); RWPxy = (time_meas,hgt);
        #######################################################################
        #RWPdata = (WS,WD,Uvel,Vvel,Wvel,SNRbeam,SWbeam,Wqual)
        RWPdata = (WS,WD,Uvel,Vvel,SNRbeam,Wqual)
        #######################################################################
        #RWPlow = (HGTlow,WSlow,WDlow,Uvel_low,Vvel_low,Wvel_low,SNRlow,SNRlow_rngecorr,SWlow,Wlow_qual)
        #RWPhgh = (HGThgh,WShgh,WDhgh,Uvel_hgh,Vvel_hgh,Wvel_hgh,SNRhgh,SNRhgh_rngecorr,SWhgh,Whgh_qual)
        RWPlow = (HGTlow,WSlow,WDlow,Uvel_low,Vvel_low,SNRlow,SNRlow_rngecorr,Wlow_qual)
        RWPhgh = (HGThgh,WShgh,WDhgh,Uvel_hgh,Vvel_hgh,SNRhgh,SNRhgh_rngecorr,Whgh_qual)

    return RWPloc, RWPxy, RWPdata, RWPlow, RWPhgh, lev2ref
