#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The Purpose of this Script is to Extract from the Various Initialization Files 
(i.e. nc) Relevant Model Data for Comparison to Observations

Script Does Not Limit the Length of the Forecast Hours -- This is Done in the 
Individual Plotting Scripts
-- In an initial iteration of this script, model forecast output outwards of 
-- +18 UTC was only considered
###############################################################################
Created on Fri Dec 13 15:07:20 2019
@author: jduncan

###############################################################################
Input MGref Depicts the Model Analysis Product that is to be Analyzed:
    Value of 0: Indicates -- Nearest Model Grid Point Was Used
    Value of 1: Indicates -- Bilinear Interpolation Was Used

<class 'netCDF4._netCDF4.Dataset'>
root group (NETCDF3_CLASSIC data model, file format NETCDF3):
    version: $Id: extract_columns.ncl,v 1.16 2018/12/05 14:47:59 Dave.Turner Exp $
    date_file_created: File created on Sun Sep 29 23:03:18 UTC 2019
    dimensions(sizes): forecast_hour(23), sites(73), height_bin(61), depth_bin(3), maxstring(61)
    
Variables:
    Information Currently Being Extracted/Stored:
    float32 lat, float32 lon, char dsite(maxstring), 
    ###########################################################################
    float32 mlat, float32 mlon, float32 distance, float32 elev, 
    ###########################################################################
    float32 forecast(forecast_hour)
    ###########################################################################
    Precipitation Information
    float32 prate0(forecast_hour,sites), float32 prate1(forecast_hour,sites):: Precipitation Rate   
    
"""

import os
import netCDF4 as netcdf
import numpy as np
############################
#from uv_to_met import uv2wd
############################
# Define Forecast Hours of Interest -- Currently Considering Those w/in XX Hrs
#FCSThr = 18; 
# forecast hours no longer limited

def pslmdl_nc(nc_dir,nc_fil):
    """ Define the netCDF File to be Processed """
    NCoi = os.path.join(nc_dir,nc_fil);
    ###########################################################################
    """ Initiate Model Analysis """
    PSLfil = netcdf.Dataset(NCoi);
    
    """ Site Information """
    PSLlon = np.asarray(PSLfil.variables['lon']) #Desired site longitude (E)
    PSLlat = np.asarray(PSLfil.variables['lat']) #Desired site latitude (N)
    #Measurement Site Name Information
    PSLref = np.asarray(PSLfil.variables['station_name']);
            
    """ Model Gridpoint Information """
    Mlon = np.asarray(PSLfil.variables['mlon']) #Closest model grid point longitude (E)
    Mlat = np.asarray(PSLfil.variables['mlat']) #Closest model grid point latitude (N)
    Melev = np.asarray(PSLfil.variables['malt']);
    #Distance Between Desired and Model Grid Points
    #Mdist = np.asarray(PSLfil.variables['distance']);
    Mdist = np.asarray(None);
    
    """ Model Data Information """
    Mmeas_hgt = np.asarray(PSLfil.variables['malt']) #m AGL
    Mfcst_tme = np.asarray(PSLfil.variables['forecast']) #Forecast Hour 
    ###########################################################################
    #HRoi = Mfcst_tme <= FCSThr;
    # Appropriately Reduce the Model Forecast Time
    #Mfcst_tme = Mfcst_tme[HRoi]
    
    #######################################
    # Precipitation Information
    #######################################
    apcp_o = np.asarray(PSLfil.variables['precipitation']) # Precipitation rate at the closest model grid location
    #Prate_int = np.asarray(PSLfil.variables['prate1']) #Precipitation rate calculated by bilinear interpolation

    """ Create Tuples to Return to Main Script """
    ###########################################################################
    # PSL Reference Information
    PSLxy = (Mfcst_tme,Mmeas_hgt)
    PSLloc = (PSLref,PSLlon,PSLlat,Mlon,Mlat,Melev,Mdist)
    ###########################################################################
    #-> Atmospheric Variable Information Returned in Sepearte Arrays for those 
    #-> with Data Relating to the Nearest Gridpoint and those with Data Derived 
    #-> from Bilinear Interpolation
    ###########################################################################
    # PSL Precipitation
    PSLpr_o = (apcp_o)
    #PSLpr_int = (Tint,Tsfc_int,Pint,Psfc_int,MXrat_int,MXrat_sfc_int,Prate_int)
    
    ###########################################################################
    PSLpr = (PSLpr_o)
    ###########################################################################
    return PSLxy, PSLloc, PSLpr
