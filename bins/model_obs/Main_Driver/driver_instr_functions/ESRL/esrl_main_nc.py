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
    float32 dlat(sites), float32 dlon(sites), |S1 dsite(sites,maxstring), 
    ###########################################################################
    float32 mlat(sites), float32 mlon(sites), float32 distance(sites), float32 elev(sites), 
    ###########################################################################
    float32 height(height_bin), float32 forecast(forecast_hour)
    ###########################################################################
    float32 p0(forecast_hour,sites,height_bin), float32 p1(forecast_hour,sites,height_bin):: Pressure
    float32 psfc0(forecast_hour,sites), float32 psfc1(forecast_hour,sites):: Pressure at the surface
    ###########################################################################
    float32 t0(forecast_hour,sites,height_bin), float32 t1(forecast_hour,sites,height_bin):: Temperature
    float32 tsfc0(forecast_hour,sites), float32 tsfc1(forecast_hour,sites): Temperature at 2-m
    ###########################################################################
    float32 rhsfc0(forecast_hour,sites), float32 rhsfc1(forecast_hour,sites):: Relative humidity at 2-m
    float32 r0(forecast_hour,sites,height_bin), float32 r1(forecast_hour,sites,height_bin):: Water vapor mixing ratio
    float32 rsfc0(forecast_hour,sites), float32 rsfc1(forecast_hour,sites):: Water vapor mixing ratio at 2-m
    ###########################################################################
    float32 wspd0(forecast_hour,sites,height_bin), float32 wspd1(forecast_hour,sites,height_bin):: Wind speed
    ###########################################################################
    float32 u0(forecast_hour,sites,height_bin), float32 u1(forecast_hour,sites,height_bin):: U-Component of Wind
    float32 v0(forecast_hour,sites,height_bin), float32 v1(forecast_hour,sites,height_bin):: V-Component of Wind
    float32 w0(forecast_hour,sites,height_bin), float32 w1(forecast_hour,sites,height_bin):: W-Component of Wind
    ###########################################################################
    float32 usfc0(forecast_hour,sites), float32 usfc1(forecast_hour,sites):: U-component of wind at 10-m
    float32 vsfc0(forecast_hour,sites), float32 vsfc1(forecast_hour,sites):: V-component of wind at 10-m
    ###########################################################################
    float32 lflux0(forecast_hour,sites), float32 lflux1(forecast_hour,sites):: Latent heat flux
    float32 sflux0(forecast_hour,sites), float32 sflux1(forecast_hour,sites):: Sensible heat flux
    float32 ustar0(forecast_hour,sites), float32 ustar1(forecast_hour,sites):: Friction velocity
    ###########################################################################
    float32 dswsfc0(forecast_hour,sites), float32 dswsfc1(forecast_hour,sites):: Downwelling SW flux at the surface
    float32 uswsfc0(forecast_hour,sites), float32 uswsfc1(forecast_hour,sites):: Upwelling SW flux at the surface
    float32 uswtoa0(forecast_hour,sites), float32 uswtoa1(forecast_hour,sites):: Upwelling SW flux at the top of atmosphere
    
    float32 dlwsfc0(forecast_hour,sites), float32 dlwsfc1(forecast_hour,sites):: Downwelling LW flux at the surface
    float32 ulwsfc0(forecast_hour,sites), float32 ulwsfc1(forecast_hour,sites):: Upwelling LW flux at the surface
    float32 ulwtoa0(forecast_hour,sites), float32 ulwtoa1(forecast_hour,sites):: Upwelling LW flux at the top of the atmosphere
    ###########################################################################
    float32 dswdbeam0(forecast_hour,sites), float32 dswdbeam1(forecast_hour,sites):: SWdown direct beam
    float32 dswdiffuse0(forecast_hour,sites), float32 dswdiffuse1(forecast_hour,sites):: SWdown diffuse
    ###########################################################################
    Clear Sky Information
    float32 dswdbeamCS0(forecast_hour,sites), float32 dswdbeamCS1(forecast_hour,sites):: SWdown direct beam assuming clear sky
    float32 dswdiffuseCS0(forecast_hour,sites), float32 dswdiffuseCS1(forecast_hour,sites):: SW diffuse assuming clear sky
    float32 dswsfcCS0(forecast_hour,sites), float32 dswsfcCS1(forecast_hour,sites):: Downwelling SW flux at surface assuming clear sky
    float32 uswsfcCS0(forecast_hour,sites), float32 uswsfcCS1(forecast_hour,sites):: Upwelling SW flux at surface assuming clear sky
    float32 ulwsfcCS0(forecast_hour,sites), float32 ulwsfcCS1(forecast_hour,sites):: Upwelling LW flux at surface assuming clear sky
    float32 dlwsfcCS0(forecast_hour,sites), float32 dlwsfcCS1(forecast_hour,sites):: Downwelling LW flux at surface assuming clear sky
    ###########################################################################
    Precipitation Information
    float32 prate0(forecast_hour,sites), float32 prate1(forecast_hour,sites):: Precipitation Rate   
    ###########################################################################
    Misc Information
    float32 cbh0(forecast_hour,sites), float32 cbh1(forecast_hour,sites)::  Cloud Base Height
    ###########################################################################
    Information Not Currently Being Extracted/Stored:
    ###########################################################################
    ###########################################################################
    float32 depth(depth_bin)
    -----> Depth (i.e. m BGL)
    float32 qv0(forecast_hour,sites,height_bin), float32 qv1(forecast_hour,sites,height_bin),
    -----> Specific Humidity 
    float32 sigma_qv(forecast_hour,sites,height_bin)
    -----> Specific Humidity Uncertainty    
    float32 qvsfc0(forecast_hour,sites), float32 qvsfc1(forecast_hour,sites)
    -----> Specific humidity at 2-m
    ###########################################################################
    float32 tv0(forecast_hour,sites,height_bin), float32 tv1(forecast_hour,sites,height_bin) 
    -----> Virtual temperature
    float32 lwc0(forecast_hour,sites,height_bin), float32 lwc1(forecast_hour,sites,height_bin)
    -----> Liquid water mixing ratio
    float32 iwp0(forecast_hour,sites,height_bin), float32 iwp1(forecast_hour,sites,height_bin)
    -----> Cloud ice
    float32 mxlen0(forecast_hour,sites,height_bin), float32 mxlen1(forecast_hour,sites,height_bin) 
    -----> Mixing length scale
    float32 totcc0(forecast_hour,sites,height_bin), float32 totcc1(forecast_hour,sites,height_bin)
    -----> Total cloud cover 
    ###########################################################################
    float32 tke0(forecast_hour,sites,height_bin), float32 tke1(forecast_hour,sites,height_bin)
    -----> Turbulent Kinetic Energy    
    float32 soilt0(forecast_hour,sites,depth_bin), float32 soilt1(forecast_hour,sites,depth_bin)
    -----> Soil temperature 
    float32 soilm0(forecast_hour,sites,depth_bin), float32 soilm1(forecast_hour,sites,depth_bin)
    -----> Volumetric soil moisture
    float32 pwv0(forecast_hour,sites), float32 pwv1(forecast_hour,sites)
    -----> Precipitable water vapor
    float32 pbl0(forecast_hour,sites), float32 pbl1(forecast_hour,sites)
    -----> Pressure at cloud base
    float32 stemp0(forecast_hour,sites), float32 stemp1(forecast_hour,sites)
    -----> Surface skin temperature
    float32 sigma_t(forecast_hour,sites,height_bin)
    -----> Temperature uncertainty  (not in new profiles)
    ###########################################################################
    float32 gflux0(forecast_hour,sites), float32 gflux1(forecast_hour,sites)
    -----> Ground heat flux
    ###########################################################################
    float32 tcc0(forecast_hour,sites), float32 tcc1(forecast_hour,sites)
    -----> Total cloud cover fraction
    float32 lcc0(forecast_hour,sites), float32 lcc1(forecast_hour,sites)
    -----> Low cloud cover fraction
    float32 mcc0(forecast_hour,sites), float32 mcc1(forecast_hour,sites)
    -----> Medium cloud cover fraction
    float32 hcc0(forecast_hour,sites), float32 hcc1(forecast_hour,sites)
    -----> High cloud cover fraction
    float32 vis0(forecast_hour,sites), float32 vis1(forecast_hour,sites)
    -----> Visibility        
    float32 crefl0(forecast_hour,sites), float32 crefl1(forecast_hour,sites)
    -----> Composite reflectivity    
    float32 lwp0(forecast_hour,sites), float32 lwp1(forecast_hour,sites)
    -----> Vertically integrated liquid (LWP)
    float32 ri0(forecast_hour,sites), float32 ri1(forecast_hour,sites)
    -----> Richardson Number    
    float32 terr0(forecast_hour,sites), float32 terr1(forecast_hour,sites)
    -----> Terrain Height    
    ###########################################################################
    float32 csnow0(forecast_hour,sites), float32 csnow1(forecast_hour,sites)
    -----> Categorical snow
    float32 snowfrac0(forecast_hour,sites), float32 snowfrac1(forecast_hour,sites)
    -----> Snow cover
    ###########################################################################
    float32 convcld0(forecast_hour,sites), float32 convcld1(forecast_hour,sites)
    -----> Convective Cloud Cover
    float32 cldamount0(forecast_hour,sites), float32 cldamount1(forecast_hour,sites)
    -----> Maximum Cloud Amount
    float32 vegtype0(forecast_hour,sites), float32 vegtype1(forecast_hour,sites)
    -----> Vegetation Type
    float32 aod0(forecast_hour,sites), float32 aod1(forecast_hour,sites)
    -----> Aerosol Optical Depth
    float32 assa0(forecast_hour,sites), float32 assa1(forecast_hour,sites)
    -----> Aerosol Single Scatter Albedo
    float32 asym0(forecast_hour,sites), float32 asym1(forecast_hour,sites)
    -----> Aerosol Asymmetry Parameter
    ###########################################################################
    float32 dpsfc0(forecast_hour,sites), float32 dpsfc1(forecast_hour,sites)
    -----> Dew Point Temperature at 2-m
    ###########################################################################
    float32 u80m0(forecast_hour,sites), float32 u80m1(forecast_hour,sites), 
    float32 v80m0(forecast_hour,sites), float32 v80m1(forecast_hour,sites)
    -----> U and V wind components of wind at 80-m
"""

import os
import netCDF4 as netcdf
import numpy as np
############################
from uv_to_met import uv2wd
from uv_to_met import uv2ws
############################
from pslutils import read_var_fill_missing
############################
# Define Forecast Hours of Interest -- Currently Considering Those w/in XX Hrs
#FCSThr = 18; 
# forecast hours no longer limited

def esrl_nc(nc_dir,nc_fil):
    """ Define the netCDF File to be Processed """
    NCoi = os.path.join(nc_dir,nc_fil);
    ###########################################################################
    """ Initiate Model Analysis """
    ESRLfil = netcdf.Dataset(NCoi);
    ESRLvars = ESRLfil.variables
    
    """ Site Information """
    ESRLlon = np.asarray(ESRLfil.variables['dlon']) #Desired site longitude (E)
    ESRLlat = np.asarray(ESRLfil.variables['dlat']) #Desired site latitude (N)
    #Measurement Site Name Information
    ESRLref = np.asarray(ESRLfil.variables['dsite']);
    ESRLref = ESRLref.astype(str); #encoding issue correction
    ESRLsite = [''.join(row).strip(',') for row in ESRLref]; del ESRLref;
            
    """ Model Gridpoint Information """
    Mlon = np.asarray(ESRLfil.variables['mlon']) #Closest model grid point longitude (E)
    Mlat = np.asarray(ESRLfil.variables['mlat']) #Closest model grid point latitude (N)
    if 'elev' in ESRLvars:
        Melev = np.asarray(ESRLfil.variables['elev']);
    else:
        Melev = Mlon*0
    #Distance Between Desired and Model Grid Points
    Mdist = np.asarray(ESRLfil.variables['distance']);
    
    """ Model Data Information """
    Mmeas_hgt = np.asarray(ESRLfil.variables['height']) #m AGL
    Mfcst_tme = np.asarray(ESRLfil.variables['forecast']) #Forecast Hour 
    ###########################################################################
    #HRoi = Mfcst_tme <= FCSThr;
    # Appropriately Reduce the Model Forecast Time
    #Mfcst_tme = Mfcst_tme[HRoi]
    
    """ Temperature, Pressure, and Relative Humidity Atmospheric Variables """ 
    ######################################################
    #Atmospheric Pressure (Units: hPa)
    Po = np.asarray(ESRLfil.variables['p0']) #Pressure at the closest model grid location
    Pint = np.asarray(ESRLfil.variables['p1']) #Pressure calculated by bilienar interpolation
    ###########################################################################
    # Reduce to Forecast Times of Interest
    #Po = Po[HRoi,:,:];
    #Pint = Pint[HRoi,:,:];
    
    Psfc_o = np.asarray(ESRLfil.variables['psfc0']) #Pressure at the surface at the closest model grid location
    Psfc_int = np.asarray(ESRLfil.variables['psfc1']) #Pressure at the surface calculated by bilinear interpolation
    ###########################################################################
    # Reduce to Forecast Times of Interest
    #Psfc_o = Psfc_o[HRoi,:]
    #Psfc_int = Psfc_int[HRoi,:];

    ###########################################################################
    #Atmospheric Temperature (Units: C)
    To = np.asarray(ESRLfil.variables['t0']) #Temperature at the closest model grid location
    Tint = np.asarray(ESRLfil.variables['t1']) #Temperature calculated by bilienar interpolation
    #new profiles don't have sigma_t
    #To_sigma = np.asarray(ESRLfil.variables['sigma_t']) #Temperature uncertainty around the closest model grid location
    #To_sigma = np.zeros_like(To) 
    ###########################################################################
    # Reduce to Forecast Times of Interest
    #To = To[HRoi,:,:];
    #Tint = Tint[HRoi,:,:];
    #To_sigma = To_sigma[HRoi,:,:];
    
    Tsfc_o = np.asarray(ESRLfil.variables['tsfc0']) #Temperature at 2-m at the closest model grid location
    Tsfc_int = np.asarray(ESRLfil.variables['tsfc1']) #Temperature at 2-m calculated by bilinear interpolation
    ###########################################################################
    
    """ Relative Humidity Information No Longer Used for Model Evaluation """
    ###########################################################################
    # Relative Humidity (Units: %)
    RHsfc_o = read_var_fill_missing(ESRLfil,'rhsfc0',Tsfc_o) #Relative humidity at 2-m at the closest model grid location
    RHsfc_int = read_var_fill_missing(ESRLfil,'rhsfc1',Tsfc_int) #Relative humidity at 2-m calculated by bilinear interpolation
    #if 'rhsfc0' in ESRLvars:
    #    RHsfc_o = np.asarray(ESRLfil.variables['rhsfc0']) #Relative humidity at 2-m at the closest model grid location
    #    RHsfc_int = np.asarray(ESRLfil.variables['rhsfc1']) #Relative humidity at 2-m calculated by bilinear interpolation
    #else:
    #    RHsfc_o = Tsfc_o*np.nan
    #    RHsfc_int = Tsfc_int*np.nan
    ###########################################################################
    
    ###########################################################################
    """ Mixing Ratio will Replace Relative Humidity for Model Evaluation """
    ###########################################################################
    # Water Vapor Mixing Ratio (Units: g/kg)
    MXrat_o = np.asarray(ESRLfil.variables['r0']) #Water vapor mixing ratio at the closest model grid location
    MXrat_int = np.asarray(ESRLfil.variables['r1']) #Water vapor mixing ratio calculated by bilinear interpolation
    # Fill Value: 9.96921e+36
    ###########################################################################
    # Surface Mixing Ratio Values (Units g/kg)
    MXrat_sfc_o = np.asarray(ESRLfil.variables['rsfc0']) #Water vapor mixing ratio at 2-m at the closest model grid location
    MXrat_sfc_int = np.asarray(ESRLfil.variables['rsfc1']) #Water vapor mixing ratio at 2-m calculated by bilinear interpolation
    # Fill Value: 9.96921e+36
    ###########################################################################

    #######################################
    # Precipitation Information (Units kg/m2s = mm/s)
    #######################################
    Prate_o = np.asarray(ESRLfil.variables['prate0']) # Precipitation rate at the closest model grid location
    Prate_int = np.asarray(ESRLfil.variables['prate1']) #Precipitation rate calculated by bilinear interpolation

    """ Accumulated Precip - added in v2.0.0 """
    ###########################################################################
    # Hourly Accumulated Precipitation (Units: kg/m2)
    if 'paccum0' in ESRLvars:
        Paccum_o = np.asarray(ESRLfil.variables['paccum0']) #Hourly accumulated precip at the closest model grid location
        Paccum_int = np.asarray(ESRLfil.variables['paccum1']) #Hourly accumulated precip calculated by bilinear interpolation
    else:
        Paccum_o = Prate_o*np.nan
        Paccum_int = Prate_int*np.nan

    #######################################
    # Precipitable Water Information (Units cm)
    #######################################
    Pwtr_o = np.asarray(ESRLfil.variables['pwv0']) # Precipitable water at the closest model grid location
    Pwtr_int = np.asarray(ESRLfil.variables['pwv1']) #Precipitable water calculated by bilinear interpolation

    #######################################
    # Planetary Boundary Layer Information (Units km)
    #######################################
    Pbl_o = np.asarray(ESRLfil.variables['pbl0']) # PBL at the closest model grid location
    Pbl_int = np.asarray(ESRLfil.variables['pbl1']) #PBL calculated by bilinear interpolation

    #######################################
    # Cloud Base Height (Units  hPa - pressure level)
    #######################################
    Cbh_o = read_var_fill_missing(ESRLfil,'cbh0',Tsfc_o) # CBH at the closest model grid location
    Cbh_int = read_var_fill_missing(ESRLfil,'cbh1',Tsfc_int) # CBH calculated by bilinear interpolation
    #if 'cbh0' in ESRLvars:
    #    Cbh_o = np.asarray(ESRLfil.variables['cbh0']) # CBH at the closest model grid location
    #    Cbh_int = np.asarray(ESRLfil.variables['cbh1']) #CBH calculated by bilinear interpolation
    #else:
    #    Cbh_o = Pbl_o*np.nan
    #    Cbh_int = Pbl_int*np.nan

    ###########################################################################
    """ Wind Speed and Direction Atmospheric Variables """
    ###########################################################################
    WSo = np.asarray(ESRLfil.variables['wspd0']) #Wind speed at the closest model grid location
    WSint = np.asarray(ESRLfil.variables['wspd1']) #Wind speed calculated by bilinear interpolation
    # Units m/s
    ###########################################################################
    
    ###########################################################################
    #U,V, and W Wind Components (Units m/s)
    Uo = np.asarray(ESRLfil.variables['u0']) #u-component of wind at the closest model grid location
    Uint = np.asarray(ESRLfil.variables['u1']) #u-component of wind calculated by bilinear interpolation
    ###########################################################################
    Vo = np.asarray(ESRLfil.variables['v0']) #v-component of wind at the closest model grid location
    Vint = np.asarray(ESRLfil.variables['v1']) #v-component of wind calculated by bilinear interpolation
    ###########################################################################
    Wo = np.asarray(ESRLfil.variables['w0']) #w-component of wind at the closest model grid location
    Wint = np.asarray(ESRLfil.variables['w1']) #w-component of wind calculated by bilinear interpolation
    ###########################################################################
    # Determine Wind Direction (Already Reduced to Proper Dimensions)
    WDo = uv2wd(Uo,Vo); WDint = uv2wd(Uint,Vint);
    #############################################
    
    # Define Surface Information
    ###########################################################################
    Usfc_o = np.asarray(ESRLfil.variables['usfc0']) #U-component of wind at 10-m at the closest model grid location
    Usfc_int = np.asarray(ESRLfil.variables['usfc1']) #U-component of wind at 10-m calculated by bilinear interpolation
    ###########################################################################
    Vsfc_o = np.asarray(ESRLfil.variables['vsfc0']) #V-component of wind at 10-m at the closest model grid location
    Vsfc_int = np.asarray(ESRLfil.variables['vsfc1']) #V-component of wind at 10-m calculated by bilinear interpolation
    ###########################################################################

    # Surface Wind Speed and Direction (Already Reduced to Proper Dimensions)
    #WSsfc_o = np.sqrt(Usfc_o**2 + Vsfc_o**2); WSsfc_int = np.sqrt(Usfc_int**2 + Vsfc_int**2);
    WSsfc_o = uv2ws(Usfc_o,Vsfc_o); WSsfc_int = uv2ws(Usfc_int,Vsfc_int);
    WDsfc_o = uv2wd(Usfc_o,Vsfc_o); WDsfc_int = uv2wd(Usfc_int,Vsfc_int);
    
    """ Flux Atmospheric Variables """
    ##################################
    # Latent Heat Flux (Units W/m2)
    LHo = np.asarray(ESRLfil.variables['lflux0']) #Latent heat flux at the closest model grid location
    LHint = np.asarray(ESRLfil.variables['lflux1']) #Latent heat flux calculated by bilinear interpolation
    ###########################################################################
    
    # Sensible Heat Flux (Units W/m2)
    So = np.asarray(ESRLfil.variables['sflux0']) #Sensible heat flux at the closest model grid location
    Sint = np.asarray(ESRLfil.variables['sflux1']) #Sensible heat flux calculated by bilinear interpolation
    ###########################################################################
    
    # Friction Velocity (Units m/s)
    Ustar_o = np.asarray(ESRLfil.variables['ustar0']) #Friction velocity at the closest model grid location
    Ustar_int = np.asarray(ESRLfil.variables['ustar1']) #Friction velocity calculated by bilinear interpolation

    """ Radiation Atmospheric Variables """
    ##########################################
    # Radiation Flux Information (Units: W/m2)
    ##########################################
    SWdsfc_o = read_var_fill_missing(ESRLfil,'dswsfc0',Tsfc_o) # Downwelling SW flux at the surface at the closest model grid location
    SWdsfc_int = read_var_fill_missing(ESRLfil,'dswsfc1',Tsfc_int) # Downwelling SW flux at the surface calculated by bilinear interpolation
    #SWdsfc_o = np.asarray(ESRLfil.variables['dswsfc0']) #Downwelling SW flux at the surface at the closest model grid location
    #SWdsfc_o[np.where(SWdsfc_o == ESRLfil.variables['dswsfc0']._FillValue)] = np.nan
    #SWdsfc_int = np.asarray(ESRLfil.variables['dswsfc1']) #Downwelling SW flux at the surface calculated by bilinear interpolation
    #SWdsfc_int[np.where(SWdsfc_int == ESRLfil.variables['dswsfc1']._FillValue)] = np.nan
    ###########################################################################
    SWusfc_o = read_var_fill_missing(ESRLfil,'uswsfc0',Tsfc_o) # Downwelling SW flux at the surface at the closest model grid location
    SWusfc_int = read_var_fill_missing(ESRLfil,'uswsfc1',Tsfc_int) # Downwelling SW flux at the surface calculated by bilinear interpolation
    #SWusfc_o = np.asarray(ESRLfil.variables['uswsfc0']) #Upwelling SW flux at the surface at the closest model grid location
    #SWusfc_o[np.where(SWusfc_o == ESRLfil.variables['uswsfc0']._FillValue)] = np.nan
    #SWusfc_int = np.asarray(ESRLfil.variables['uswsfc1']) #Upwelling SW flux at the surface calculated by bilinear interpolation
    #SWusfc_int[np.where(SWusfc_int == ESRLfil.variables['uswsfc1']._FillValue)] = np.nan
    ###########################################################################
    SWutoa_o = read_var_fill_missing(ESRLfil,'uswtoa0',Tsfc_o) # Downwelling SW flux at the surface at the closest model grid location
    SWutoa_int = read_var_fill_missing(ESRLfil,'uswtoa1',Tsfc_int) # Downwelling SW flux at the surface calculated by bilinear interpolation
    #SWutoa_o = np.asarray(ESRLfil.variables['uswtoa0']) #Upwelling SW flux at the top of atmosphere at the closest model grid location
    #SWutoa_o[np.where(SWutoa_o == ESRLfil.variables['uswtoa0']._FillValue)] = np.nan
    #SWutoa_int = np.asarray(ESRLfil.variables['uswtoa1']) #Upwelling SW flux at the top of atmosphere calculated by bilinear interpolation
    #SWutoa_int[np.where(SWutoa_int == ESRLfil.variables['uswtoa1']._FillValue)] = np.nan
    ###########################################################################
    LWdsfc_o = np.asarray(ESRLfil.variables['dlwsfc0']) #Downwelling LW flux at the surface at the closest model grid location
    LWdsfc_o[np.where(LWdsfc_o == ESRLfil.variables['dlwsfc0']._FillValue)] = np.nan
    LWdsfc_int = np.asarray(ESRLfil.variables['dlwsfc1']) #Downwelling LW flux at the surface calculated by bilinear interpolation
    LWdsfc_int[np.where(LWdsfc_int == ESRLfil.variables['dlwsfc1']._FillValue)] = np.nan
    ###########################################################################
    LWusfc_o = np.asarray(ESRLfil.variables['ulwsfc0']) #Upwelling LW flux at the surface at the closest model grid location
    LWusfc_o[np.where(LWusfc_o == ESRLfil.variables['ulwsfc0']._FillValue)] = np.nan
    LWusfc_int = np.asarray(ESRLfil.variables['ulwsfc1']) #Upwelling LW flux at the surface calculated by bilinear interpolation
    LWusfc_int[np.where(LWusfc_int == ESRLfil.variables['ulwsfc1']._FillValue)] = np.nan
    ###########################################################################
    LWutoa_o = np.asarray(ESRLfil.variables['ulwtoa0']) #Upwelling LW flux at the top of atmosphere at the closest model grid location
    LWutoa_o[np.where(LWutoa_o == ESRLfil.variables['ulwtoa0']._FillValue)] = np.nan
    LWutoa_int = np.asarray(ESRLfil.variables['ulwtoa1']) #Upwelling LW flux at the top of atmosphere calculated by bilinear interpolation
    LWutoa_int[np.where(LWutoa_int == ESRLfil.variables['ulwtoa1']._FillValue)] = np.nan
    ###########################################################################
    
    #######################################
    # Direct Radiation Information
    #######################################
    if 'dswdbeam0' in ESRLvars:
        SWddir_o = np.asarray(ESRLfil.variables['dswdbeam0']) #SWdown direct beam at the closest model grid location
        SWddir_int = np.asarray(ESRLfil.variables['dswdbeam1']) #SWdown direct beam calculated by bilinear interpolation
        SWddif_o = np.asarray(ESRLfil.variables['dswdiffuse0']) #SWdown diffuse at the closest model grid location
        SWddif_int = np.asarray(ESRLfil.variables['dswdiffuse1']) #SWdown diffuse calculated by bilinear interpolation
    else:
        SWddir_o = SWdsfc_o*np.nan
        SWddir_int = SWdsfc_int*np.nan
        SWddif_o = SWdsfc_o*np.nan
        SWddif_int = SWdsfc_int*np.nan
    ###########################################################################
    
    #######################################
    # Clear Sky Radiation Information
    #######################################
# not used and not in the new GSL profiles
### SWddircs_o = np.asarray(ESRLfil.variables['dswdbeamCS0']) #SWdown direct beam assuming clear sky at the closest model grid location
### SWddircs_int = np.asarray(ESRLfil.variables['dswdbeamCS1']) #SWdown direct beam assuming clear sky calculated by bilinear interpolation
### 
### SWddifcs_o = np.asarray(ESRLfil.variables['dswdiffuseCS0']) #SW diffuse assuming clear sky at the closest model grid location
### SWddifcs_int = np.asarray(ESRLfil.variables['dswdiffuseCS1']) #SW diffuse assuming clear sky calculated by bilinear interpolation
### 
### SWdsfccs_o = np.asarray(ESRLfil.variables['dswsfcCS0']) #Downwelling SW flux at surface assuming clear sky at the closest model grid location
### SWdsfccs_int = np.asarray(ESRLfil.variables['dswsfcCS1']) #Downwelling SW flux at surface assuming clear sky calculated by bilinear interpolation
### 
### SWusfccs_o = np.asarray(ESRLfil.variables['uswsfcCS0']) #Upwelling SW flux at surface assuming clear sky at the closest model grid location
### SWusfccs_int = np.asarray(ESRLfil.variables['uswsfcCS1']) #Upwelling SW flux at surface assuming clear sky calculated by bilinear interpolation
### 
### LWusfccs_o = np.asarray(ESRLfil.variables['ulwsfcCS0']) #Upwelling LW flux at surface assuming clear sky at the closest model grid location
### LWusfccs_int = np.asarray(ESRLfil.variables['ulwsfcCS1']) #Upwelling LW flux at surface assuming clear sky calculated by bilinear interpolation
### 
### LWdsfccs_o = np.asarray(ESRLfil.variables['dlwsfcCS0']) #Downwelling LW flux at surface assuming clear sky at the closest model grid location
### LWdsfccs_int = np.asarray(ESRLfil.variables['dlwsfcCS1']) #Downwelling LW flux at surface assuming clear sky calculated by bilinear interpolation
### ###########################################################################
### # Reduce to Forecast Times of Interest
### #SWddircs_o = SWddircs_o[HRoi,:]; SWddircs_int = SWddircs_int[HRoi,:]; 
### #SWddifcs_o = SWddifcs_o[HRoi,:]; SWddifcs_int = SWddifcs_int[HRoi,:]; 
### #SWdsfccs_o = SWdsfccs_o[HRoi,:]; SWdsfccs_int = SWdsfccs_int[HRoi,:]; 
### 
### #SWusfccs_o = SWusfccs_o[HRoi,:]; SWusfccs_int = SWusfccs_int[HRoi,:]; 
### #LWusfccs_o = LWusfccs_o[HRoi,:]; LWusfccs_int = LWusfccs_int[HRoi,:]; 
### #LWdsfccs_o = LWdsfccs_o[HRoi,:]; LWdsfccs_int = LWdsfccs_int[HRoi,:]; 

    """ Create Tuples to Return to Main Script """
    ###########################################################################
    # ESRL Reference Information
    ESRLxy = (Mfcst_tme,Mmeas_hgt)
    ESRLloc = (ESRLsite,ESRLlon,ESRLlat,Mlon,Mlat,Melev,Mdist)
    ###########################################################################
    #-> Atmospheric Variable Information Returned in Sepearte Arrays for those 
    #-> with Data Relating to the Nearest Gridpoint and those with Data Derived 
    #-> from Bilinear Interpolation
    ###########################################################################
    # ESRL Temperature, Pressure, Mixing Ratio Information, Precipitation, PWV, PBL
    ESRLtpmx_o = (To,Tsfc_o,Po,Psfc_o,MXrat_o,MXrat_sfc_o,Prate_o,Pwtr_o,Pbl_o,Cbh_o,RHsfc_o,Paccum_o)
    ESRLtpmx_int = (Tint,Tsfc_int,Pint,Psfc_int,MXrat_int,MXrat_sfc_int,Prate_int,Pwtr_int,Pbl_int,Cbh_int,RHsfc_int,Paccum_int)
    
    ###########################################################################
    ESRLtpmx = (ESRLtpmx_o,ESRLtpmx_int)
    ###########################################################################
    # ESRL Wind Speed, Direction, and Wind Components
    ESRLw_o = (WSo,WDo,WSsfc_o,WDsfc_o,Uo,Usfc_o,Vo,Vsfc_o,Wo)
    ESRLw_int = (WSint,WDint,WSsfc_int,WDsfc_int,Uint,Usfc_int,Vint,Vsfc_int,Wint)   
    ###########################################################################
    ESRLw = (ESRLw_o,ESRLw_int)
    ###########################################################################
    # ESRL Flux Information            
    ESRLflux_o = (LHo,So,Ustar_o)
    ESRLflux_int = (LHint, Sint, Ustar_int)
    ###########################################################################
    ESRLflux = (ESRLflux_o,ESRLflux_int)
    ###########################################################################
    # ESRL Radiation Information
    ESRLradfluxsfc_o = (SWdsfc_o, LWdsfc_o,SWusfc_o, LWusfc_o)
    ESRLradfluxsfc_int = (SWdsfc_int, LWdsfc_int, SWusfc_int, LWusfc_int)
    
    ESRLradfluxtoa_o = (SWutoa_o,LWutoa_o)
    ESRLradfluxtoa_int = (SWutoa_int,LWutoa_int)
    ###########################################################################
    ESRLrad_flux = (ESRLradfluxsfc_o,ESRLradfluxtoa_o,ESRLradfluxsfc_int,ESRLradfluxtoa_int)
    ###########################################################################
    ESRLdd_o = (SWddir_o,SWddif_o)
    ESRLdd_int = (SWddir_int,SWddif_int)
    ###########################################################################
    ESRLrad_dd = (ESRLdd_o, ESRLdd_int)    
    ###########################################################################
    #ESRLradcs_o = (SWdsfccs_o,LWdsfccs_o,SWusfccs_o,LWusfccs_o,SWddircs_o,SWddifcs_o)
    #ESRLradcs_int = (SWdsfccs_int,LWdsfccs_int,SWusfccs_int,LWusfccs_int,SWddircs_int,SWddifcs_int)  
    ###########################################################################
    #ESRLrad_cs = (ESRLradcs_o,ESRLradcs_int)
    ESRLrad_cs = ()
    ###########################################################################
    ESRLrad  = (ESRLrad_flux, ESRLrad_dd, ESRLrad_cs)
    ###########################################################################
    return ESRLxy, ESRLloc, ESRLtpmx, ESRLw, ESRLflux, ESRLrad
