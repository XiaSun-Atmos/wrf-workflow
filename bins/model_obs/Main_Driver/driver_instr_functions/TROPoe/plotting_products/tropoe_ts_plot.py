#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This Function Plots for the Case where both Measurement and Model Data are 
Available

###############################################################################
Created on Wed Jan 31 2023
@author: dmurray
###############################################################################
The Purpose of this Script is to Develop and Output time series plots
Plots Based on TROPoe retrievals
###############################################################################
Outputted Plots Follow:
instrument_whatisplotted_Sid.date.png
**Imperative that the "." is the first period -- this is used to strip the **
**date from the file so that duplicate images are not run **
###############################################################################
Most Plots Generated in This Script Will Consist of Two Panels 
- in the first plot, the obs will be plotted along with the model output
- in the second plot, the model bias will be shown and statistics will be annotated

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as image
import os,time,warnings
from datetime import datetime,timedelta
#####################################
from coord_distance import dist2coord
from atmo_meas_to_hour_end import atmo_meas2hr_end
from time_tuple_cat import time_cat
from atmo_tuple_cat import atmo_cat
from atmo_time_correction import atmo_time_correct
##########################################
from atmo_bias import atmo_bias
from atmo_ylim_bias import atmo_blim
from atmo_ylim_plot import atmo_ylim
from atmo_ylim_bounded import atmo_bound_lim
from label_funcs import avgmins_to_label
from label_funcs import gridpoint_to_label
from format_funcs import format_lat_lon_alt
####################################
#Initiate Plotting Parameter Module
import plotting_parameters as plt_param
import station_param_list as station_param
import forecast_param_list as fxcst_param

################################################################
# Define Model Interpolation Paramter
INTPoi = plt_param.INToi;
#0 indicates nearest-neighbor
#1 indicates bilinear interpolation
###############################################################################
HRsec = 60*60;
DAYsec = 24*HRsec;
# averaging period for comparison (minutes), label for plot
AVGper = getattr(plt_param,'TROPoe_avgmins',60)
AVGlabel = avgmins_to_label(AVGper)
DEFpblh=[300,318]
###############################################################################

def tropoe_tsplot(MEASin,MDLin,OUTinfo,version='ASSIST'):
    DATE = OUTinfo[3];
    basedate=datetime.strptime(DATE,'%Y%m%d')
    ###########################################################################
    IMGpth = os.path.join(OUTinfo[0],DATE);
    if not os.path.exists(IMGpth):
        os.makedirs(IMGpth)
    IMGpth_ext = os.path.join(OUTinfo[0],'ext',DATE);
    if not os.path.exists(IMGpth_ext):
        os.makedirs(IMGpth_ext)
    IMGpth_model = os.path.join(OUTinfo[0],'no_model',DATE);
    if not os.path.exists(IMGpth_model):
        os.makedirs(IMGpth_model)
    ###########################################################################Z
    MDLname = OUTinfo[1];
    Sid = OUTinfo[2];
    logo=image.imread(plt_param.NOAAlogo);
    FILpfx='tropoe'+version.lower()+'_'
    if 'ASSIST' in version:
       ATTRstr = "TROPoe Retrieval Infrared Spectrometer"
    else:
       ATTRstr = "TROPoe Retrieval Microwave Radiometer"
    
    """ Based on MDLname -- Define Relevant Variables """
    FXlen = getattr(fxcst_param,MDLname.lower()+'_fx')
    EXTlen = getattr(fxcst_param,MDLname.lower()+'_ext')
    #if MDLname == 'HRRR_v4': # 'ESRL_HRRR';
    #    FXlen = fxcst_param.hrrr_v4_fx;
    #    EXTlen = fxcst_param.hrrr_v4_ext;
    #elif MDLname == 'RAP_v5': # 'ESRL_RAP';
    #    FXlen = fxcst_param.rap_v5_fx;
    #    EXTlen = fxcst_param.rap_v5_ext;
    # +1 Hr to Forecast Length to Account for Ini Z Forecast Inclusion
    ###########################################################################
    """ Explicit Definition of the Input Tuples
    tropoe_TSplot(TROPoe_plt,HRRRw_plt,FILout)
    #############################################
    OUTinfo :: FILout = (IMGout,MDLname,Sid,doi)
    #############################################"""
    """ Unpack Modeled T/P/MX Data
    ###########################
    MDLin :: 
            MDLin = (MDL_ini,MDL_loc,MDL_xy,MDL_tpmx), e.g., 
              HRRRin = (HRRRv4_ini,HRRRv4_loc,HRRRv4_xy,HRRRv4_tpmx)
              RAPin = (RAPv5_ini,RAPv5_loc,RAPv5_xy,RAPv5_tpmx)
            #-> each embedded tuple has length of initialization files
            ###################################################################
            MDLloc = (PSLsite,PSLlon,PSLlat,Mlon,Mlat,Melev,Mdist)
            MDLxy = (Mfcst_tme,Mmeas_hgt)
            MDLtpmx = depends on interpolation choice (see below)
                    ESRLtpmx_o = (To,Tsfc_o,Po,Psfc_o,MXrat_o,MXrat_sfc_o,Prate_o,Pwv_o,Pblh_o)
                    ESRLtpmx_int = (Tint,Tsfc_int,Pint,Psfc_int,MXrat_int,MXrat_sfc_int,Prate_int,Pwv_int,Pblh_int)
                    ###########################################################################
                    ESRLtpmx = (ESRLtpmx_o,ESRLtpmx_int)
    ####################################################################### """
    # Model Initialization Times #        
    MDLini = MDLin[0]; #outputs model initialization time for each ini hour
    MDLxy = MDLin[2]; #outputs forecast time and measurement height
    MDLtpmx = MDLin[3] #outputs tpmx ini hour (expand in loop)
    ###########################################################################
    MDLloc = MDLin[1][0]; #just store first initialization hour
    MDLloc_names = MDLloc[0];
    MDLlon_sites = MDLloc[1]; MDLlat_sites = MDLloc[2]; # desired grid point
    ######################################################################
    MDLlon_grids = MDLloc[3]; MDLlat_grids = MDLloc[4];# model grid point referenced
    MDLelev_grids = MDLloc[5]
    MDLdist2sites = MDLloc[6]
    ###########################################################################

    ###########################################################################
    """ Define Measurement Location and Model Grid Point Structure """
    MEASlon = getattr(station_param,Sid+'lon')
    MEASlat = getattr(station_param,Sid+'lat')
    MEASalt = getattr(station_param,Sid+'alt')
    MEASid = getattr(station_param,Sid+'name')
    MEASname = getattr(station_param,Sid+'descr')
    ###########################################################################
    # Get alternate grid point if available
    MEASmdlid = getattr(station_param,Sid+'mdlID',"") # mdl ID for multiple points
    ALTsite = getattr(station_param,Sid+'altSite',None)
    ALTmdlid = getattr(station_param,Sid+'altmdlID',"") # mdl ID for alternate point
    ALTmdls = getattr(station_param,Sid+'altmdls',[]) # mdl ID for alternate point
    #MEAStitle = MEASname+' ('+str(MEASalt)+' m)'
    MEAStitle = MEASname+'\n'+format_lat_lon_alt(MEASlat,MEASlon,MEASalt)
    ###########################################################################
    """ Unpack TROPoe Information 
    ####################################
    MEASin :: TROPoe_plt = (TROPoeloc,TROPoetme,TROPoetrh,TROPoetheta,TROPoemisc)
    ---> TROPoeloc = (TROPoelon,TROPoelat,TROPoealt)
    ---> TROPoetme = (TROPoexy_minus,TROPoexy,TROPoexy_plus,TROPoexy_ext)
    -----> TROPoetrh = (TROPoe_trh_minus,TROPoe_trh,TROPoe_trh_plus,TROPoe_trh_ext)
    -----> TROPoetheta = (TROPoe_theta_minus,TROPoe_theta,TROPoe_theta_plus,TROPoe_theta_ext)
    -----> TROPoemisc = (TROPoe_misc_minus,TROPoe_misc,TROPoe_misc_plus,TROPoe_misc_ext)
    ###########################################################################
    TROPoeloc = (TROPoelon,TROPoelat,TROPoealt); TROPoexy = (time_meas,hgt);
    ###########################################################################
    TROPoe_trh = (temp,rh,dewpt,mx,pres)
    TROPoe_theta = (theta,thetae)
    TROPoe_misc = (cbh,gamma,rmsr,lwp,pwv,pblh)
    ########################################################################"""
    # Define Relevant Indices for Reference in Concatenation #
    CBHind_cat = 0; GAMMAind_cat = 1; RMSRind_cat=2; LWPind_cat=3; PWVind_cat=4; PBLind_cat=5;
    
    #-> Combine Relevant Measurement Tuples 
    MEAStime = time_cat(MEASin[1][0],MEASin[1][1],MEASin[1][2],MEASin[1][3])
    ###########################################################################
    MEASmisc_in = MEASin[4]; 
    ###########################################################################
    #MEAScbh = atmo_cat(MEASmisc_in[0],MEASmisc_in[1],MEASmisc_in[2],MEASmisc_in[3],CBHind_cat)*1000 # to meters
    MEASgamma = atmo_cat(MEASmisc_in[0],MEASmisc_in[1],MEASmisc_in[2],MEASmisc_in[3],GAMMAind_cat)
    MEASrmsr = atmo_cat(MEASmisc_in[0],MEASmisc_in[1],MEASmisc_in[2],MEASmisc_in[3],RMSRind_cat)
    MEASlwp = atmo_cat(MEASmisc_in[0],MEASmisc_in[1],MEASmisc_in[2],MEASmisc_in[3],LWPind_cat)
    MEASpwv = atmo_cat(MEASmisc_in[0],MEASmisc_in[1],MEASmisc_in[2],MEASmisc_in[3],PWVind_cat)*10 # to mm
    MEASpblh = atmo_cat(MEASmisc_in[0],MEASmisc_in[1],MEASmisc_in[2],MEASmisc_in[3],PBLind_cat)*1000 # to meters
    
    # set values to nan when certain conditions are met
    for t in range(len(MEASgamma)):
       if not np.logical_and(MEASgamma[t] == 1, MEASrmsr[t] < 5):
           MEASpwv[t] = np.nan
           MEASpblh[t] = np.nan
           MEASlwp[t] = np.nan

    # Per Bianca Adler, LWP > 200 is unrealistic so set it to NaN
    MEASlwp[np.where(MEASlwp > 200)] = np.nan

    # don't use pblh when the value is close to the default values
    for dval in DEFpblh:
        MEASpblh[np.where(np.isclose(MEASpblh,dval,atol=.5))] = np.nan;

    #  For ASSIST plots, only use data when lwp <= 5
    if (version == 'ASSIST'):
        MEASpwv[np.where(MEASlwp > 5)] = np.nan;
        MEASpblh[np.where(MEASlwp > 5)] = np.nan;

    ###################################################################
    # Perform Correction to Create Equally Spaced Data Arrays #
    ###################################################################
    MEAStime_new, MEASpwv = atmo_time_correct(MEAStime,MEASpwv);
    _, MEASpblh = atmo_time_correct(MEAStime,MEASpblh);
    _, MEASlwp = atmo_time_correct(MEAStime,MEASlwp);
    MEAStime = MEAStime_new

    # Identify the Nearest Model Grid Point to the Site Of Interest
    Grefs = []
    if (INTPoi == 0):
        Gref_meas, Gref_dist = dist2coord((MEASlat,MEASlon),(MDLlat_grids,MDLlon_grids))
    else:
        Gref_meas, Gref_dist = dist2coord((MEASlat,MEASlon),(MDLlat_sites,MDLlon_sites))
    Grefs.append(Gref_meas)
    ###########################################################################
    # If an alternate site has been requested, find the model data closest to that.
    if (ALTsite and (MDLname in ALTmdls)):
       ALTlat = ALTsite[0]; ALTlon = ALTsite[1];
       if (INTPoi == 0):
           Gref_alt, Gref_dist = dist2coord((ALTlat,ALTlon),(MDLlat_grids,MDLlon_grids))
       else:
           Gref_alt, Gref_dist = dist2coord((ALTlat,ALTlon),(MDLlat_sites,MDLlon_sites))
       Grefs.append(Gref_alt)
    ###########################################################################
    # Set up some labelling for the plots based on the sites selected
    sidSfx = "";   # suffix for the station ID in the plot name
    sfxNum = 0;     # the suffix counter
    MDLstr = ""     # qualifier for the model location for plot labelling
    for Gref in Grefs:
        if (sfxNum > 0):   # the MEAS station has no suffix
           sidSfx = str(sfxNum)
        Sid = Sid+sidSfx
        #  add in the qualifier based on which point we are processing
        if (ALTsite and (MDLname in ALTmdls)): 
            if (sfxNum==0):
                MDLstr = ' ('+MEASmdlid+')'
            else:
                MDLstr = ' ('+ALTmdlid+')'

        ###########################################################################
        # Store the Reference Grid Information
        MDLlon_site =  MDLlon_sites[Gref]; MDLlat_site = MDLlat_sites[Gref]; 
        MDLlon_grid =  MDLlon_grids[Gref]; MDLlat_grid = MDLlat_grids[Gref];   
        MDLelev_grid = MDLelev_grids[Gref]
        MDLdist2site = MDLdist2sites[Gref];
        MDLloc_name = MDLloc_names[Gref]
        MDLgridinfo = gridpoint_to_label(INTPoi,sfxNum,Grefs,MDLlat_grid,MDLlon_grid,MDLelev_grid,\
                                         MDLlat_site,MDLlon_site,MEASmdlid,ALTmdlid)
        ###########################################################################
        sfxNum = sfxNum+1
            
        # Define the Range of Possible Initialization Times #
        INIpos = ['%02d' % it for it in range(0,24)]
        """ Initiate Measurement Model Comparison for Each Initialization Hour """
        for ini in range(0,len(INIpos)):
            # Determine Whether Model Data Exists for Said Date
            if INIpos[ini] in MDLini:
                """ Extract Relevant Model Data """
                ###########################################
                INIref = MDLini.index(INIpos[ini])
                ###########################################
                # Unpack Relevant Model Information
                MDLtime = MDLxy[INIref][0];
                ###################################################################
                # -> Convert Time (Initialization Time + Forecast Hour) to Seconds 
                # -> After Midnight
                ###################################################################
                MDLtime = int(MDLini[INIref])*HRsec + MDLtime*HRsec
                ###################################################################
    
                """ Extract Relevant Model Information """
                MDLtpmx_ini = MDLtpmx[INIref][INTPoi]
                ###################################################################
                # Extract Estimates of precipitable water and boundary layer height
                MDLpwv = MDLtpmx_ini[7]
                MDLpblh = MDLtpmx_ini[8]
                # surface--not height-values were extracted
                ###################################################################
                # Extract from the Relevant Grid Point
                MDLpwv = MDLpwv[:,Gref];
                MDLpblh = MDLpblh[:,Gref];
                # until we have real LWP values, set them all to NaN
                MDLlwp = MDLpblh*np.nan
                ###################################################################
                """ Define Relevant Plotting Specifics 
                --> Specifics will vary depending on whether an extended forecast 
                --> is being produced
                ###################################################################
                Extended is Necessary if Forecast Exceeds Model-Specific Standard 
                Length
                """
                if len(MDLtime) > FXlen: #indicating forecast hours exceed standard forecast length
                    """ Define Bounds for Extended Forecast Period """
                    ext_meas_oi = np.logical_and(MEAStime >= MDLtime[0], MEAStime <= MDLtime[0] + EXTlen*HRsec) 
                    ###############################################################
                    # Adjust Model Variables #
                    MDLpwv_ext = MDLpwv; MDLpblh_ext = MDLpblh; MDLlwp_ext = MDLlwp;
                    MDLtime_ext = MDLtime;
                    ###############################################################
                    # Perform Model Averaging for Plotting and Statistical Analsyes
                    MEASpwv_ext_mu = atmo_meas2hr_end(MEASpwv,None,MEAStime,MDLtime_ext,AVGmins=AVGper);
                    MEASpblh_ext_mu = atmo_meas2hr_end(MEASpblh,None,MEAStime,MDLtime_ext,AVGmins=AVGper);
                    MEASlwp_ext_mu = atmo_meas2hr_end(MEASlwp,None,MEAStime,MDLtime_ext,AVGmins=AVGper);
                    ###############################################################
                    #PWVmin_ext,PWVmax_ext,PWVtick_ext = atmo_ylim((MEASpwv[ext_meas_oi],\
                    #                               MDLpwv_ext),plt_param.PWVrnd_bse,plt_param.TICKint)
                    #PBLHmin_ext,PBLHmax_ext,PBLHtick_ext = atmo_ylim((MEASpblh[ext_meas_oi],\
                    #                               MDLpblh_ext),plt_param.PBLHrnd_bse,plt_param.TICKint)
                    PWVmin_ext,PWVmax_ext,PWVtick_ext = atmo_bound_lim((MEASpwv[ext_meas_oi],\
                                                   MDLpwv_ext),plt_param.PWVrnd_bse,plt_param.PWVmin,plt_param.PWVmax,\
                                                   plt_param.PWVrnd_bse)
                    PBLHmin_ext,PBLHmax_ext,PBLHtick_ext = atmo_bound_lim((MEASpblh[ext_meas_oi],\
                                                   MDLpblh_ext),plt_param.PBLHrnd_bse,plt_param.PBLHmin,plt_param.PBLHmax,\
                                                   plt_param.PBLHrnd_bse)
                    LWPmin_ext,LWPmax_ext,LWPtick_ext = atmo_bound_lim((MEASlwp[ext_meas_oi],\
                                                   MDLlwp_ext),plt_param.LWPrnd_bse,plt_param.LWPmin,plt_param.LWPmax,\
                                                   plt_param.LWPrnd_bse)
                    ###############################################################
                    #Define Xtick Information (Push Out to Maximum Forecast Time)
                    Xmin_ext = MDLtime_ext[0] - HRsec;
                    Xmax_ext = (MDLtime_ext[0]+EXTlen*HRsec) + HRsec;
                    # do not use indice approach because we plan to push out fxcst 
                    ###############################################################
                    Xtick_ext = np.arange(MDLtime_ext[0],MDLtime_ext[0]+EXTlen*HRsec+1e-5,plt_param.Xdel_ext)
                    Xtick_ext_lab = [time.strftime('%H', time.gmtime(time_ref)) for time_ref in Xtick_ext]
                    Time_range_dates_ext = [basedate+timedelta(seconds=(time_ref)) for time_ref in Xtick_ext]
                    Date_start_label_ext=Time_range_dates_ext[0].strftime('%Y-%m-%d')
                    Date_end_label_ext=Time_range_dates_ext[-1].strftime('%Y-%m-%d')
    
                    """ Define Bounds for the Model-Specific Standard Forecast Period"""
                    meas_oi = np.logical_and(MEAStime >= MDLtime[0], MEAStime <= MDLtime[FXlen-1]) 
                    ###############################################################
                    # Adjust Model Variables #
                    MDLpwv = MDLpwv[0:FXlen]; MDLpblh = MDLpblh[0:FXlen]; MDLlwp = MDLlwp[0:FXlen];
                    MDLtime = MDLtime[0:FXlen];
                    ###############################################################
                    # Perform Model Averaging for Plotting and Statistical Analsyes
                    MEASpwv_mu = atmo_meas2hr_end(MEASpwv,None,MEAStime,MDLtime,AVGmins=AVGper);
                    MEASpblh_mu = atmo_meas2hr_end(MEASpblh,None,MEAStime,MDLtime,AVGmins=AVGper);
                    MEASlwp_mu = atmo_meas2hr_end(MEASlwp,None,MEAStime,MDLtime,AVGmins=AVGper);
                    ###############################################################
                    #PWVmin,PWVmax,PWVtick = atmo_ylim((MEASpwv[meas_oi],\
                    #                               MDLpwv),plt_param.PWVrnd_bse,plt_param.TICKint)
                    #PBLHmin,PBLHmax,PBLHtick = atmo_ylim((MEASpblh[meas_oi],\
                    #                               MDLpblh),plt_param.PBLHrnd_bse,plt_param.TICKint)
                    PWVmin,PWVmax,PWVtick = atmo_bound_lim((MEASpwv[meas_oi],\
                                                   MDLpwv),plt_param.PWVrnd_bse,plt_param.PWVmin,plt_param.PWVmax,\
                                                   plt_param.PWVrnd_bse)
                    PBLHmin,PBLHmax,PBLHtick = atmo_bound_lim((MEASpblh[meas_oi],\
                                                   MDLpblh),plt_param.PBLHrnd_bse,plt_param.PBLHmin,plt_param.PBLHmax,\
                                                   plt_param.PBLHrnd_bse)
                    LWPmin,LWPmax,LWPtick = atmo_bound_lim((MEASlwp[meas_oi],\
                                                   MDLlwp),plt_param.LWPrnd_bse,plt_param.LWPmin,plt_param.LWPmax,\
                                                   plt_param.LWPrnd_bse)
                    ###############################################################
                    #Define Xtick Information (Push Out to Maximum Forecast Time)
                    Xmin = MDLtime[0] - HRsec;
                    Xmax = (MDLtime[0]+(FXlen - 1)*HRsec) + HRsec;
                    # do not use indice approach because we plan to push out fxcst 
                    ###############################################################
                    Xtick = np.arange(MDLtime[0],MDLtime[0]+(FXlen - 1)*HRsec+1e-5,plt_param.Xdel)
                    Xtick_lab = [time.strftime('%H', time.gmtime(time_ref)) for time_ref in Xtick]
                    Time_range_dates = [basedate+timedelta(seconds=(time_ref)) for time_ref in Xtick]
                    Date_start_label=Time_range_dates[0].strftime('%Y-%m-%d')
                    Date_end_label=Time_range_dates[-1].strftime('%Y-%m-%d')
                    
                else:
                    """ Define Bounds for the Model-Specific Standard Forecast Period"""
                    ###############################################################
                    # Know that an Extended Forecast is Non-Existent. Therefore,
                    # Simply Define meas_oi based on the Desired Forecast Time of 
                    # 86400
                    ###############################################################
                    meas_oi = np.logical_and(MEAStime >= MDLtime[0], MEAStime <= MDLtime[0]+(FXlen - 1)*HRsec)
                    ###############################################################
                    # no model adjustment needed -- either equal to + (FXlen - 1) Fxcst or
                    # or less due to a download issue (i.e. MDLws = MDLws[0:FXlen])
                    ###############################################################
                    # Perform Model Averaging for Plotting and Statistical Analsyes
                    MEASpwv_mu = atmo_meas2hr_end(MEASpwv,None,MEAStime,MDLtime,AVGmins=AVGper);
                    MEASpblh_mu = atmo_meas2hr_end(MEASpblh,None,MEAStime,MDLtime,AVGmins=AVGper);
                    MEASlwp_mu = atmo_meas2hr_end(MEASlwp,None,MEAStime,MDLtime,AVGmins=AVGper);
                    ###############################################################
                    #PWVmin,PWVmax,PWVtick = atmo_ylim((MEASpwv[meas_oi],\
                    #                               MDLpwv),plt_param.PWVrnd_bse,plt_param.TICKint)
                    #PBLHmin,PBLHmax,PBLHtick = atmo_ylim((MEASpblh[meas_oi],\
                    #                               MDLpblh),plt_param.PBLHrnd_bse,plt_param.TICKint)
                    PWVmin,PWVmax,PWVtick = atmo_bound_lim((MEASpwv[meas_oi],\
                                                   MDLpwv),plt_param.PWVrnd_bse,plt_param.PWVmin,plt_param.PWVmax,\
                                                   plt_param.PWVrnd_bse)
                    PBLHmin,PBLHmax,PBLHtick = atmo_bound_lim((MEASpblh[meas_oi],\
                                                   MDLpblh),plt_param.PBLHrnd_bse,plt_param.PBLHmin,plt_param.PBLHmax,\
                                                   plt_param.PBLHrnd_bse)
                    LWPmin,LWPmax,LWPtick = atmo_bound_lim((MEASlwp[meas_oi],\
                                                   MDLlwp),plt_param.LWPrnd_bse,plt_param.LWPmin,plt_param.LWPmax,\
                                                   plt_param.LWPrnd_bse)
                    ###############################################################
                    #Define Xtick Information (Push Out to Maximum Forecast Time)
                    Xmin = MDLtime[0] - HRsec;
                    Xmax = (MDLtime[0]+(FXlen - 1)*HRsec) + HRsec;
                    # do not use indice approach because we plan to push out fxcst 
                    ###############################################################
                    Xtick = np.arange(MDLtime[0],MDLtime[0]+(FXlen - 1)*HRsec+1e-5,plt_param.Xdel)
                    Xtick_lab = [time.strftime('%H', time.gmtime(time_ref)) for time_ref in Xtick]
                    Time_range_dates = [basedate+timedelta(seconds=(time_ref)) for time_ref in Xtick]
                    Date_start_label=Time_range_dates[0].strftime('%Y-%m-%d')
                    Date_end_label=Time_range_dates[-1].strftime('%Y-%m-%d')
    
                ###################################################################
                """ Measurements and Model Data are Available """
                """ ###############################################################
                Precipitable Water Profile (instruct to share x-axis)
                ############################################################### """
                pwv_fig, (pwv_ax,bias_ax) = plt.subplots(2,1,sharex = True, figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                         facecolor = 'w',edgecolor ='w',clear='True')
                ####################################################################
                # Plot the Precipitable Water Time Histories
                #####################################
                #Observations -- Only Plotting the QC'd Values
                pwv_ax.plot(MDLtime,MEASpwv_mu,'s',color='gray',markersize=plt_param.PLTmsize,label='Obs ('+AVGlabel+' Avg)')
                pwv_ax.plot(MEAStime,MEASpwv, color='black',linewidth=plt_param.PLTlwidth,label='Obs')
                pwv_ax.plot(MEAStime,MEASpwv,'.',color='black',markersize=1,linestyle='None',label='Obs')
                pwv_ax.plot(MDLtime,MDLpwv,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname+MDLstr)
                pwv_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ####################################################################
                handles, labels = pwv_ax.get_legend_handles_labels()
                order = [1,0,3]
                pwv_ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], \
                               loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                               ncol = 3,frameon = False,columnspacing = plt_param.LEGcol_space);    
                ###################################################################
                pwv_ax.set_ylim(PWVmin,PWVmax); pwv_ax.set_yticks(PWVtick);
                pwv_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                pwv_ax.set_ylabel('Precipitable Water \nVapor, mm',{'fontsize': plt_param.LABELfs});
                pwv_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                ###################################################################
                #insert station info as title
                T = pwv_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                T.set_position(plt_param.Tloc);
                             
                """ ###############################################################
                # Plot the Model Precipitable Water Biases
                ############################################################### """
                # initate 0/0 line
                bias_ax.plot([Xmin,Xmax],[0,0],color='black',linewidth=plt_param.PLTlwidth)
                ###################################################################
                # Define the Bias (Observations - Measurements)
                PWVbias = MDLpwv - MEASpwv_mu
                bias_ax.plot(MDLtime,PWVbias,'^-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Mdl - Obs')
                bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                bias_ax.legend(loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                             ncol = 1,frameon = False);    
                ###################################################################
                bias_ax.set_xticks(Xtick); bias_ax.set_xticklabels(Xtick_lab,fontsize=plt_param.TICKfs);
                bias_ax.set_xlim(Xmin,Xmax)
                ###################################################################
                bias_ax.set_xlabel('Hour, UTC',{'fontsize': plt_param.LABELfs});
                bias_ax.xaxis.set_label_coords(plt_param.Xlab_pad_dual[0],plt_param.Xlab_pad_dual[1]);
                ####################################################################
                #Define Bounds of Precipitable Water - Model Bias to Determine Y-Bounds
                #PWVb_min,PWVb_max,PWVb_tick = atmo_blim((PWVbias),plt_param.PWVbias_bse,plt_param.TICKint)
                PWVb_min,PWVb_max,PWVb_tick = atmo_bound_lim((PWVbias),plt_param.PWVbias_bse,plt_param.PWVbmin,\
                                                     plt_param.PWVbmax,plt_param.PWVbias_bse)
                bias_ax.set_ylim(PWVb_min,PWVb_max); bias_ax.set_yticks(PWVb_tick);
                bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                bias_ax.set_ylabel('Model Error, mm',{'fontsize': plt_param.LABELfs});    
                bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                ###################################################################
                #input basic statistics
                RMSEpwv,BIASpwv = atmo_bias(PWVbias)
                #insert data statistics as title
                BT = bias_ax.set_title('RMSE: ' + str(RMSEpwv)+' mm | Mean Bias: '+\
                                      str(BIASpwv)+' mm',{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5) 
                BT.set_position(plt_param.Bloc_left);
                ###################################################################
                pwv_fig.subplots_adjust(left = 0.125,right = 0.935,bottom = 0.155, top = 0.925)
                ###################################################################
                bias_ax.annotate(Date_start_label,xy = plt_param.metDATEpad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                bias_ax.annotate(Date_end_label,xy = plt_param.metDATE_end_pad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                #######################################################################
                # Add some extra info
                bias_ax.annotate(MDLgridinfo,xy = plt_param.GRIDpad_left_dual,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'left',\
                                fontsize = plt_param.TITLEfs)
                bias_ax.annotate(ATTRstr,xy = plt_param.ATTRpad_ctr_dual,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'center',\
                                fontsize = plt_param.TITLEfs)
                #######################################################################
                # add a logo
                imax = pwv_fig.add_axes(plt_param.metLOGOaxcoords)
                imax.set_axis_off()
                imax.imshow(logo, aspect="equal")
                #######################################################################
                #Save the Developed Figure
                FILout = os.path.join(IMGpth,FILpfx+MDLname+'_PWV_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                pwv_fig.savefig(FILout,dpi=300); plt.close('all');
    
                """ ###############################################################
                #Planetary Boundary Layer Plots (instruct to share x-axis)
                ############################################################### """
                pblh_fig, (pblh_ax,bias_ax) = plt.subplots(2,1,sharex = True, figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                         facecolor = 'w',edgecolor ='w',clear='True')
                ###################################################################
                # Plot the Atmospheric Planetary Boundary Layer Time Histories
                ##############################################
                #Observations 
                pblh_ax.plot(MDLtime,MEASpblh_mu,'s',color='gray',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Obs ('+AVGlabel+' Avg)')
                pblh_ax.plot(MEAStime,MEASpblh, color='black',linewidth=plt_param.PLTlwidth,label='Obs')
                pblh_ax.plot(MEAStime,MEASpblh,'.',color='black',markersize=1,linestyle='None',label='Obs')
                pblh_ax.plot(MDLtime,MDLpblh,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname+MDLstr)
                pblh_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                handles, labels = pblh_ax.get_legend_handles_labels()
                order = [1,0,3]
                pblh_ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], \
                               loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                               ncol = 3,frameon = False,columnspacing = plt_param.LEGcol_space);    
                ###################################################################
                pblh_ax.set_ylim(PBLHmin,PBLHmax); pblh_ax.set_yticks(PBLHtick);
                pblh_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                pblh_ax.set_ylabel('Boundary Layer Height,\n m AGL',{'fontsize': plt_param.LABELfs});
                pblh_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                ###################################################################
                #insert station info as title
                T = pblh_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                T.set_position(plt_param.Tloc);
                
                """ ###############################################################
                # Plot the Model Boundary Layer Biases
                ##############################################################  """
                # initate 0/0 line
                bias_ax.plot([Xmin,Xmax],[0,0],color='black',linewidth=plt_param.PLTlwidth)
                ###################################################################
                # Define the Bias (Observations - Measurements)
                PBLHbias = MDLpblh - MEASpblh_mu
                bias_ax.plot(MDLtime,PBLHbias,'^-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Mdl - Obs')
                bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                bias_ax.legend(loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                             ncol = 1,frameon = False);    
                ###################################################################
                bias_ax.set_xticks(Xtick); bias_ax.set_xticklabels(Xtick_lab,fontsize=plt_param.TICKfs);
                bias_ax.set_xlim(Xmin,Xmax)
                ###################################################################
                bias_ax.set_xlabel('Hour, UTC',{'fontsize': plt_param.LABELfs});
                bias_ax.xaxis.set_label_coords(plt_param.Xlab_pad_dual[0],plt_param.Xlab_pad_dual[1]);
                ###################################################################
                #Define Bounds of Boundary Layer - Model Bias to Determine Y-Bounds
                #PBLHb_min,PBLHb_max,PBLHb_tick = atmo_blim((PBLHbias),plt_param.PBLHbias_bse,plt_param.TICKint)
                PBLHb_min,PBLHb_max,PBLHb_tick = atmo_bound_lim((PBLHbias),plt_param.PBLHbias_bse,plt_param.PBLHbmin,\
                                                     plt_param.PBLHbmax,plt_param.PBLHbias_bse)
                bias_ax.set_ylim(PBLHb_min,PBLHb_max); bias_ax.set_yticks(PBLHb_tick);
                bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                bias_ax.set_ylabel('Model Error, m',{'fontsize': plt_param.LABELfs});    
                bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                ###################################################################
                #input basic statistics
                RMSEpblh,BIASpblh = atmo_bias(PBLHbias)
                #insert data statistics as title
                BT = bias_ax.set_title('RMSE: ' + str(RMSEpblh)+' m | Mean Bias: '+\
                                      str(BIASpblh)+' m',{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5) 
                BT.set_position(plt_param.Bloc_left);
                ###################################################################
                pblh_fig.subplots_adjust(left = 0.125,right = 0.935,bottom = 0.155, top = 0.925)
                ###################################################################
                bias_ax.annotate(Date_start_label,xy = plt_param.metDATEpad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                bias_ax.annotate(Date_end_label,xy = plt_param.metDATE_end_pad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                #######################################################################
                # Add some extra info
                bias_ax.annotate(MDLgridinfo,xy = plt_param.GRIDpad_left_dual,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'left',\
                                fontsize = plt_param.TITLEfs)
                bias_ax.annotate(ATTRstr,xy = plt_param.ATTRpad_ctr_dual,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'center',\
                                fontsize = plt_param.TITLEfs)
                #######################################################################
                # add a logo
                imax = pblh_fig.add_axes(plt_param.metLOGOaxcoords)
                imax.set_axis_off()
                imax.imshow(logo, aspect="equal")
                #######################################################################
                #Save the Developed Figure
                FILout = os.path.join(IMGpth,FILpfx+MDLname+'_PBLH_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                pblh_fig.savefig(FILout,dpi=300); plt.close('all');
                ###################################################################
                
                """ ###############################################################
                #Liquid Water Path Plots (instruct to share x-axis)
                ############################################################### """
                lwp_fig, (lwp_ax,bias_ax) = plt.subplots(2,1,sharex = True, figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                         facecolor = 'w',edgecolor ='w',clear='True')
                ###################################################################
                # Plot the Liquid Water Path Time Histories
                ##############################################
                #Observations 
                lwp_ax.plot(MDLtime,MEASlwp_mu,'s',color='gray',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Obs ('+AVGlabel+' Avg)')
                lwp_ax.plot(MEAStime,MEASlwp, color='black',linewidth=plt_param.PLTlwidth,label='Obs')
                lwp_ax.plot(MEAStime,MEASlwp,'.',color='black',markersize=1,linestyle='None',label='Obs')
                lwp_ax.plot(MDLtime,MDLlwp,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname+MDLstr)
                lwp_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                handles, labels = lwp_ax.get_legend_handles_labels()
                order = [1,0,3]
                lwp_ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], \
                               loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                               ncol = 3,frameon = False,columnspacing = plt_param.LEGcol_space);    
                ###################################################################
                lwp_ax.set_ylim(LWPmin,LWPmax); lwp_ax.set_yticks(LWPtick);
                lwp_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                lwp_ax.set_ylabel('Liquid Water Path,\n g m-2',{'fontsize': plt_param.LABELfs});
                lwp_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                ###################################################################
                #insert station info as title
                T = lwp_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                T.set_position(plt_param.Tloc);
                
                """ ###############################################################
                # Plot the Liquid Water Path Biases
                ##############################################################  """
                # initate 0/0 line
                bias_ax.plot([Xmin,Xmax],[0,0],color='black',linewidth=plt_param.PLTlwidth)
                ###################################################################
                # Define the Bias (Observations - Measurements)
                if (np.isnan(MDLlwp).all()):
                    bias_ax.plot(MDLtime,MDLtime*np.nan,'^-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Mdl - Obs')
                else:
                    LWPbias = MDLlwp - MEASlwp_mu
                    bias_ax.plot(MDLtime,LWPbias,'^-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Mdl - Obs')
                bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                bias_ax.legend(loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                             ncol = 1,frameon = False);    
                ###################################################################
                bias_ax.set_xticks(Xtick); bias_ax.set_xticklabels(Xtick_lab,fontsize=plt_param.TICKfs);
                bias_ax.set_xlim(Xmin,Xmax)
                ###################################################################
                bias_ax.set_xlabel('Hour, UTC',{'fontsize': plt_param.LABELfs});
                bias_ax.xaxis.set_label_coords(plt_param.Xlab_pad_dual[0],plt_param.Xlab_pad_dual[1]);
                ###################################################################
                #Define Bounds of Boundary Layer - Model Bias to Determine Y-Bounds
                if (np.isnan(MDLlwp).all()):
                    LWPb_min=-3; LWPb_max = 3; LWPb_tick = np.arange(LWPb_min, LWPb_max+1, 1);
                else:
                    LWPb_min,LWPb_max,LWPb_tick = atmo_bound_lim((LWPbias),plt_param.LWPbias_bse,plt_param.LWPbmin,\
                                                     plt_param.LWPbmax,plt_param.LWPbias_bse)
                bias_ax.set_ylim(LWPb_min,LWPb_max); bias_ax.set_yticks(LWPb_tick);
                bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                bias_ax.set_ylabel('Model Error, g m-2',{'fontsize': plt_param.LABELfs});    
                bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                ###################################################################
                #input basic statistics
                if (np.isnan(MDLlwp).all()):
                    # Denote that Model Data was Not Available
                    bias_ax.text(np.nanmean([Xmin,Xmax]),0,MDLname+' Data Not Available',\
                         horizontalalignment='center',verticalalignment='bottom',color='red',\
                         fontweight='semibold',fontsize=plt_param.TITLEfs);
                    BT = bias_ax.set_title('RMSE: N/A',{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                else:
                    RMSElwp,BIASlwp = atmo_bias(LWPbias)
                    #insert data statistics as title
                    BT = bias_ax.set_title('RMSE: ' + str(RMSElwp)+' m | Mean Bias: '+\
                                          str(BIASlwp)+' m',{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5) 
                BT.set_position(plt_param.Bloc_left);
                ###################################################################
                lwp_fig.subplots_adjust(left = 0.125,right = 0.935,bottom = 0.155, top = 0.925)
                ###################################################################
                bias_ax.annotate(Date_start_label,xy = plt_param.metDATEpad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                bias_ax.annotate(Date_end_label,xy = plt_param.metDATE_end_pad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                #######################################################################
                # Add some extra info
                bias_ax.annotate(MDLgridinfo,xy = plt_param.GRIDpad_left_dual,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'left',\
                                fontsize = plt_param.TITLEfs)
                bias_ax.annotate(ATTRstr,xy = plt_param.ATTRpad_ctr_dual,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'center',\
                                fontsize = plt_param.TITLEfs)
                #######################################################################
                # add a logo
                imax = lwp_fig.add_axes(plt_param.metLOGOaxcoords)
                imax.set_axis_off()
                imax.imshow(logo, aspect="equal")
                #######################################################################
                #Save the Developed Figure
                FILout = os.path.join(IMGpth,FILpfx+MDLname+'_LWP_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                lwp_fig.savefig(FILout,dpi=300); plt.close('all');
                ###################################################################
                """ Plot Extended Forecast Period if Applicable """
                ###################################################################
                if 'MDLtime_ext' in locals():
                    """ Extended Precipitable Water Forecast """
                    pwv_fig, (pwv_ax,bias_ax) = plt.subplots(2,1,sharex = True, figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                         facecolor = 'w',edgecolor ='w',clear='True')
                    ###############################################################
                    # Plot the Precipitable Water Time Histories
                    #####################################
                    #Observations 
                    pwv_ax.plot(MDLtime_ext,MEASpwv_ext_mu,'s',color='gray',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Obs ('+AVGlabel+' Avg)')
                    pwv_ax.plot(MEAStime,MEASpwv, color='black',linewidth=plt_param.PLTlwidth,label='Obs')
                    pwv_ax.plot(MEAStime,MEASpwv,'.',color='black',markersize=1,linestyle='None',label='Obs')
                    pwv_ax.plot(MDLtime_ext,MDLpwv_ext,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname+MDLstr)
                    pwv_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                    ###############################################################
                    handles, labels = pwv_ax.get_legend_handles_labels()
                    order = [1,0,3]
                    pwv_ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], \
                                   loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                                   ncol = 3,frameon = False,columnspacing = plt_param.LEGcol_space);    
                    ###############################################################
                    pwv_ax.set_ylim(PWVmin_ext,PWVmax_ext); pwv_ax.set_yticks(PWVtick_ext);
                    pwv_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                    ###############################################################
                    pwv_ax.set_ylabel('Precipitable Water \nVapor, mm',{'fontsize': plt_param.LABELfs});
                    pwv_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                    ###############################################################
                    #insert station info as title
                    T = pwv_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                    T.set_position(plt_param.Tloc);
                                 
                    """ ###########################################################
                    # Plot the Model Precipitable Water Biases
                    ########################################################### """
                    # initate 0/0 line
                    bias_ax.plot([Xmin_ext,Xmax_ext],[0,0],color='black',linewidth=plt_param.PLTlwidth)
                    ###############################################################
                    # Define the Bias (Observations - Measurements)
                    PWVbias = MDLpwv_ext - MEASpwv_ext_mu
                    bias_ax.plot(MDLtime_ext,PWVbias,'^-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Mdl - Obs')
                    bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                    ###############################################################
                    bias_ax.legend(loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                                 ncol = 1,frameon = False);    
                    ###############################################################
                    bias_ax.set_xticks(Xtick_ext); bias_ax.set_xticklabels(Xtick_ext_lab,fontsize=plt_param.TICKfs);
                    bias_ax.set_xlim(Xmin_ext,Xmax_ext)
                    ###############################################################
                    bias_ax.set_xlabel('Hour, UTC',{'fontsize': plt_param.LABELfs});
                    bias_ax.xaxis.set_label_coords(plt_param.Xlab_pad_dual[0],plt_param.Xlab_pad_dual[1]);
                    ###############################################################
                    #Define Bounds of Precipitable Water - Model Bias to Determine Y-Bounds
                    #PWVb_min,PWVb_max,PWVb_tick = atmo_blim((PWVbias),plt_param.PWVbias_bse,plt_param.TICKint)
                    PWVb_min,PWVb_max,PWVb_tick = atmo_bound_lim((PWVbias),plt_param.PWVbias_bse,plt_param.PWVbmin,\
                                                     plt_param.PWVbmax,plt_param.PWVbias_bse)
                    bias_ax.set_ylim(PWVb_min,PWVb_max); bias_ax.set_yticks(PWVb_tick);
                    bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                    ###############################################################
                    bias_ax.set_ylabel('Model Error, mm',{'fontsize': plt_param.LABELfs});    
                    bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                    ###############################################################
                    #input basic statistics
                    RMSEpwv,BIASpwv = atmo_bias(PWVbias)
                    #insert data statistics as title
                    BT = bias_ax.set_title('RMSE: ' + str(RMSEpwv)+' mm | Mean Bias: '+\
                                          str(BIASpwv)+' mm',{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5) 
                    BT.set_position(plt_param.Bloc_left);
                    ###############################################################
                    pwv_fig.subplots_adjust(left = 0.125,right = 0.935,bottom = 0.155, top = 0.925)
                    ###############################################################
                    bias_ax.annotate(Date_start_label_ext,xy = plt_param.metDATEpad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                    bias_ax.annotate(Date_end_label_ext,xy = plt_param.metDATE_end_pad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                    #######################################################################
                    # Add some extra info
                    bias_ax.annotate(MDLgridinfo,xy = plt_param.GRIDpad_left_dual,xycoords = 'figure fraction',\
                                    verticalalignment = 'bottom',horizontalalignment = 'left',\
                                    fontsize = plt_param.TITLEfs)
                    bias_ax.annotate(ATTRstr,xy = plt_param.ATTRpad_ctr_dual,xycoords = 'figure fraction',\
                                    verticalalignment = 'bottom',horizontalalignment = 'center',\
                                    fontsize = plt_param.TITLEfs)
                    #######################################################################
                    # add a logo
                    imax = pwv_fig.add_axes(plt_param.metLOGOaxcoords)
                    imax.set_axis_off()
                    imax.imshow(logo, aspect="equal")
                    #######################################################################
                    #Save the Developed Figure
                    FILout = os.path.join(IMGpth_ext,FILpfx+MDLname+'_PWV_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                    pwv_fig.savefig(FILout,dpi=300); plt.close('all');
                    ###############################################################
                    
                    """ Extended Boundary Layer Forecasts """
                    pblh_fig, (pblh_ax,bias_ax) = plt.subplots(2,1,sharex = True, figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                         facecolor = 'w',edgecolor ='w',clear='True')
                    ###############################################################
                    # Plot the Atmospheric Boundary Layer Time Histories
                    ##############################################
                    #Observations -- Only Plotting the QC'd Values
                    pblh_ax.plot(MDLtime_ext,MEASpblh_ext_mu,'s',color='gray',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Obs ('+AVGlabel+' Avg)')
                    pblh_ax.plot(MEAStime,MEASpblh, color='black',linewidth=plt_param.PLTlwidth,label='Obs')
                    pblh_ax.plot(MEAStime,MEASpblh,'.',color='black',markersize=1,linestyle='None',label='Obs')
                    pblh_ax.plot(MDLtime_ext,MDLpblh_ext,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname+MDLstr)
                    pblh_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                    ###############################################################
                    handles, labels = pblh_ax.get_legend_handles_labels()
                    order = [1,0,3]
                    pblh_ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], \
                                   loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                                   ncol = 3,frameon = False,columnspacing = plt_param.LEGcol_space);    
                    ###############################################################
                    pblh_ax.set_ylim(PBLHmin_ext,PBLHmax_ext); pblh_ax.set_yticks(PBLHtick_ext);
                    pblh_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                    ###############################################################
                    pblh_ax.set_ylabel('Boundary Layer Height, \nm AGL',{'fontsize': plt_param.LABELfs});
                    pblh_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                    ###############################################################
                    #insert station info as title
                    T = pblh_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                    T.set_position(plt_param.Tloc);
                    
                    """ ###########################################################
                    # Plot the Model Boundary Layer Biases
                    ##########################################################  """
                    # initate 0/0 line
                    bias_ax.plot([Xmin_ext,Xmax_ext],[0,0],color='black',linewidth=plt_param.PLTlwidth)
                    ###############################################################
                    # Define the Bias (Observations - Measurements)
                    PBLHbias = MDLpblh_ext - MEASpblh_ext_mu
                    bias_ax.plot(MDLtime_ext,PBLHbias,'^-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Mdl - Obs')
                    bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                    ###############################################################
                    bias_ax.legend(loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                                 ncol = 1,frameon = False);    
                    ###############################################################
                    bias_ax.set_xticks(Xtick_ext); bias_ax.set_xticklabels(Xtick_ext_lab,fontsize=plt_param.TICKfs);
                    bias_ax.set_xlim(Xmin_ext,Xmax_ext)
                    ###############################################################
                    bias_ax.set_xlabel('Hour, UTC',{'fontsize': plt_param.LABELfs});
                    bias_ax.xaxis.set_label_coords(plt_param.Xlab_pad_dual[0],plt_param.Xlab_pad_dual[1]);
                    ###############################################################
                    #Define Bounds of Boundary Layer - Model Bias to Determine Y-Bounds
                    #PBLHb_min,PBLHb_max,PBLHb_tick = atmo_blim((PBLHbias),plt_param.PBLHbias_bse,plt_param.TICKint)
                    PBLHb_min,PBLHb_max,PBLHb_tick = atmo_bound_lim((PBLHbias),plt_param.PBLHbias_bse,plt_param.PBLHbmin,\
                                                     plt_param.PBLHbmax,plt_param.PBLHbias_bse)
                    bias_ax.set_ylim(PBLHb_min,PBLHb_max); bias_ax.set_yticks(PBLHb_tick);
                    bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                    ###############################################################
                    bias_ax.set_ylabel('Model Error, m',{'fontsize': plt_param.LABELfs});    
                    bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                    ###############################################################
                    #input basic statistics
                    RMSEpblh,BIASpblh = atmo_bias(PBLHbias)
                    #insert data statistics as title
                    BT = bias_ax.set_title('RMSE: ' + str(RMSEpblh)+' m | Mean Bias: '+\
                                          str(BIASpblh)+' m',{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5) 
                    BT.set_position(plt_param.Bloc_left);
                    ###############################################################
                    pblh_fig.subplots_adjust(left = 0.125,right = 0.935,bottom = 0.155, top = 0.925)
                    ###############################################################
                    bias_ax.annotate(Date_start_label_ext,xy = plt_param.metDATEpad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                    bias_ax.annotate(Date_end_label_ext,xy = plt_param.metDATE_end_pad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                    #######################################################################
                    # Add some extra info
                    bias_ax.annotate(MDLgridinfo,xy = plt_param.GRIDpad_left_dual,xycoords = 'figure fraction',\
                                    verticalalignment = 'bottom',horizontalalignment = 'left',\
                                    fontsize = plt_param.TITLEfs)
                    bias_ax.annotate(ATTRstr,xy = plt_param.ATTRpad_ctr_dual,xycoords = 'figure fraction',\
                                    verticalalignment = 'bottom',horizontalalignment = 'center',\
                                    fontsize = plt_param.TITLEfs)
                    #######################################################################
                    # add a logo
                    imax = pblh_fig.add_axes(plt_param.metLOGOaxcoords)
                    imax.set_axis_off()
                    imax.imshow(logo, aspect="equal")
                    #######################################################################
                    #Save the Developed Figure
                    FILout = os.path.join(IMGpth_ext,FILpfx+MDLname+'_PBLH_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                    pblh_fig.savefig(FILout,dpi=300); plt.close('all');
    
                    """ Extended Liquid Water Path Forecasts """
                    lwp_fig, (lwp_ax,bias_ax) = plt.subplots(2,1,sharex = True, figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                         facecolor = 'w',edgecolor ='w',clear='True')
                    ###############################################################
                    # Plot the Liquid Water Path Time Histories
                    ##############################################
                    #Observations -- Only Plotting the QC'd Values
                    lwp_ax.plot(MDLtime_ext,MEASlwp_ext_mu,'s',color='gray',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Obs ('+AVGlabel+' Avg)')
                    lwp_ax.plot(MEAStime,MEASlwp, color='black',linewidth=plt_param.PLTlwidth,label='Obs')
                    lwp_ax.plot(MEAStime,MEASlwp,'.',color='black',markersize=1,linestyle='None',label='Obs')
                    lwp_ax.plot(MDLtime_ext,MDLlwp_ext,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname+MDLstr)
                    lwp_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                    ###############################################################
                    handles, labels = lwp_ax.get_legend_handles_labels()
                    order = [1,0,3]
                    lwp_ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], \
                                   loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                                   ncol = 3,frameon = False,columnspacing = plt_param.LEGcol_space);    
                    ###############################################################
                    lwp_ax.set_ylim(LWPmin_ext,LWPmax_ext); lwp_ax.set_yticks(LWPtick_ext);
                    lwp_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                    ###############################################################
                    lwp_ax.set_ylabel('Liquid Water Path, \ng m-2',{'fontsize': plt_param.LABELfs});
                    lwp_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                    ###############################################################
                    #insert station info as title
                    T = lwp_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                    T.set_position(plt_param.Tloc);
                    
                    """ ###########################################################
                    # Plot the Liquid Water Path Layer Biases
                    ##########################################################  """
                    # initate 0/0 line
                    bias_ax.plot([Xmin_ext,Xmax_ext],[0,0],color='black',linewidth=plt_param.PLTlwidth)
                    ###############################################################
                    # Define the Bias (Observations - Measurements)
                    if (np.isnan(MDLlwp).all()):
                        bias_ax.plot(MDLtime,MDLtime*np.nan,'^-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Mdl - Obs')
                    else:
                        LWPbias = MDLlwp_ext - MEASlwp_ext_mu
                        bias_ax.plot(MDLtime_ext,LWPbias,'^-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Mdl - Obs')
                    bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                    ###############################################################
                    bias_ax.legend(loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                                 ncol = 1,frameon = False);    
                    ###############################################################
                    bias_ax.set_xticks(Xtick_ext); bias_ax.set_xticklabels(Xtick_ext_lab,fontsize=plt_param.TICKfs);
                    bias_ax.set_xlim(Xmin_ext,Xmax_ext)
                    ###############################################################
                    bias_ax.set_xlabel('Hour, UTC',{'fontsize': plt_param.LABELfs});
                    bias_ax.xaxis.set_label_coords(plt_param.Xlab_pad_dual[0],plt_param.Xlab_pad_dual[1]);
                    ###############################################################
                    #Define Bounds of Boundary Layer - Model Bias to Determine Y-Bounds
                    if (np.isnan(MDLlwp).all()):
                        LWPb_min=-3; LWPb_max = 3; LWPb_tick = np.arange(LWPb_min, LWPb_max+1, 1);
                    else:
                        LWPb_min,LWPb_max,LWPb_tick = atmo_bound_lim((LWPbias),plt_param.LWPbias_bse,plt_param.LWPbmin,\
                                                         plt_param.LWPbmax,plt_param.LWPbias_bse)
                    bias_ax.set_ylim(LWPb_min,LWPb_max); bias_ax.set_yticks(LWPb_tick);
                    bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                    ###############################################################
                    bias_ax.set_ylabel('Model Error, g m-2',{'fontsize': plt_param.LABELfs});    
                    bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                    ###############################################################
                    #input basic statistics
                    if (np.isnan(MDLlwp).all()):
                        # Denote that Model Data was Not Available
                        bias_ax.text(np.nanmean([Xmin_ext,Xmax_ext]),0,MDLname+' Data Not Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
                        BT = bias_ax.set_title('RMSE: N/A',{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                    else:
                        RMSElwp,BIASlwp = atmo_bias(LWPbias)
                        #insert data statistics as title
                        BT = bias_ax.set_title('RMSE: ' + str(RMSElwp)+' m | Mean Bias: '+\
                                              str(BIASlwp)+' m',{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5) 
                    BT.set_position(plt_param.Bloc_left);
                    ###############################################################
                    lwp_fig.subplots_adjust(left = 0.125,right = 0.935,bottom = 0.155, top = 0.925)
                    ###############################################################
                    bias_ax.annotate(Date_start_label_ext,xy = plt_param.metDATEpad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                    bias_ax.annotate(Date_end_label_ext,xy = plt_param.metDATE_end_pad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                    #######################################################################
                    # Add some extra info
                    bias_ax.annotate(MDLgridinfo,xy = plt_param.GRIDpad_left_dual,xycoords = 'figure fraction',\
                                    verticalalignment = 'bottom',horizontalalignment = 'left',\
                                    fontsize = plt_param.TITLEfs)
                    bias_ax.annotate(ATTRstr,xy = plt_param.ATTRpad_ctr_dual,xycoords = 'figure fraction',\
                                    verticalalignment = 'bottom',horizontalalignment = 'center',\
                                    fontsize = plt_param.TITLEfs)
                    #######################################################################
                    # add a logo
                    imax = lwp_fig.add_axes(plt_param.metLOGOaxcoords)
                    imax.set_axis_off()
                    imax.imshow(logo, aspect="equal")
                    #######################################################################
                    #Save the Developed Figure
                    FILout = os.path.join(IMGpth_ext,FILpfx+MDLname+'_LWP_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                    lwp_fig.savefig(FILout,dpi=300); plt.close('all');
    
                    ###############################################################
                    # Delete Extended Information to Preclude Incorrect Plotting
                    del MDLtime_ext
                    
            else:
                """ Initialization Time Not Available """
                ###################################################################
                MDLpseudo_time = np.arange(0,(FXlen - 1)*HRsec+1e-5,step=HRsec)
                MDLtime = int(INIpos[ini])*HRsec + MDLpseudo_time
                ###################################################################
                #Define Xtick Information (Only Needs to be Defined Once)
                Xmin = MDLtime[0] - HRsec;
                Xmax = (MDLtime[0]+(FXlen - 1)*HRsec) + HRsec;
                ###################################################################
                #Define Xtick Information (Only Needs to be Defined Once)
                Xtick = np.arange(MDLtime[0],MDLtime[-1]+1e-5,plt_param.Xdel)
                Xtick_lab = [time.strftime('%H', time.gmtime(time_ref)) for time_ref in Xtick]
                Time_range_dates = [basedate+timedelta(seconds=(time_ref)) for time_ref in Xtick]
                Date_start_label=Time_range_dates[0].strftime('%Y-%m-%d')
                Date_end_label=Time_range_dates[-1].strftime('%Y-%m-%d')
                
                ###################################################################
                """ Perform Model Averaging for Plotting and Statistical Analyses """
                MEASpwv_mu = atmo_meas2hr_end(MEASpwv,None,MEAStime,MDLtime,AVGmins=AVGper);
                MEASpblh_mu = atmo_meas2hr_end(MEASpblh,None,MEAStime,MDLtime,AVGmins=AVGper);
                MEASlwp_mu = atmo_meas2hr_end(MEASlwp,None,MEAStime,MDLtime,AVGmins=AVGper);
                ###################################################################
                """ Define Relevant Plotting Specifics """
                # Define Proper Indices to Determine Axes Parameters
                meas_oi = np.logical_and(MEAStime >= MDLtime[0], MEAStime <= MDLtime[-1]) 
                ###################################################################
                #PWVmin,PWVmax,PWVtick = atmo_ylim((MEASpwv[meas_oi]),\
                #                               plt_param.PWVrnd_bse,plt_param.TICKint)
                #PBLHmin,PBLHmax,PBLHtick = atmo_ylim((MEASpblh[meas_oi]),\
                #                               plt_param.PBLHrnd_bse,plt_param.TICKint)
                PWVmin,PWVmax,PWVtick = atmo_bound_lim((MEASpwv[meas_oi]),\
                                               plt_param.PWVrnd_bse,plt_param.PWVmin,plt_param.PWVmax,\
                                               plt_param.PWVrnd_bse)
                PBLHmin,PBLHmax,PBLHtick = atmo_bound_lim((MEASpblh[meas_oi]),\
                                               plt_param.PBLHrnd_bse,plt_param.PBLHmin,plt_param.PBLHmax,\
                                               plt_param.PBLHrnd_bse)
                LWPmin,LWPmax,LWPtick = atmo_bound_lim((MEASlwp[meas_oi]),\
                                               plt_param.LWPrnd_bse,plt_param.LWPmin,plt_param.LWPmax,\
                                               plt_param.LWPrnd_bse)
                """ ###############################################################
                Precipitable Water Profile (instruct to share x-axis)
                ############################################################### """
                pwv_fig, (pwv_ax,bias_ax) = plt.subplots(2,1,sharex = True, figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                         facecolor = 'w',edgecolor ='w',clear='True')
                ####################################################################
                # Plot the Precipitable Water Time Histories
                #####################################
                #Observations 
                pwv_ax.plot(MDLtime,MEASpwv_mu,'s',color='gray',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Obs ('+AVGlabel+' Avg)')
                pwv_ax.plot(MEAStime,MEASpwv, color='black',linewidth=plt_param.PLTlwidth,label='Obs')
                pwv_ax.plot(MEAStime,MEASpwv,'.',color='black',markersize=1,linestyle='None',label='Obs')
                pwv_ax.plot(MDLtime,MDLtime*np.nan,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname+MDLstr)
                pwv_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ####################################################################
                handles, labels = pwv_ax.get_legend_handles_labels()
                order = [1,0,3]
                pwv_ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], \
                               loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                               ncol = 3,frameon = False,columnspacing = plt_param.LEGcol_space);    
                ###################################################################
                pwv_ax.set_ylim(PWVmin,PWVmax); pwv_ax.set_yticks(PWVtick);
                pwv_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                pwv_ax.set_ylabel('Precipitable Water \nVapor, mm',{'fontsize': plt_param.LABELfs});
                pwv_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                ###################################################################
                #insert data date as title
                T = pwv_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                T.set_position(plt_param.Tloc);
                             
                """ ###############################################################
                # Plot the Model Precitable Water Biases
                ############################################################### """
                # initate 0/0 line
                bias_ax.plot([Xmin,Xmax],[0,0],color='black',linewidth=plt_param.PLTlwidth)
                ###################################################################
                bias_ax.plot(MDLtime,MDLtime*np.nan,'^-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Mdl - Obs')
                bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                bias_ax.legend(loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                             ncol = 1,frameon = False);  
                ###################################################################
                #bias_ax.set_xticks(Xtick); bias_ax.set_xticklabels(Xtick_lab,{'fontsize': plt_param.TICKfs});
                bias_ax.set_xticks(Xtick); bias_ax.set_xticklabels(Xtick_lab,fontsize=plt_param.TICKfs);
                bias_ax.set_xlim(Xmin,Xmax)
                ###################################################################
                bias_ax.set_xlabel('Hour, UTC',{'fontsize': plt_param.LABELfs});
                bias_ax.xaxis.set_label_coords(plt_param.Xlab_pad_dual[0],plt_param.Xlab_pad_dual[1]);
                ###################################################################
                # Denote that Model Data was Not Availabel
                bias_ax.text(np.nanmean([Xmin,Xmax]),0.25,MDLname+' Data Not Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
                ###################################################################
                #Define Bounds of Wind Speed - Model Bias to Determine Y-Bounds
                bias_ax.set_ylim(-3,3); bias_ax.set_yticks(np.arange(-3,3.1,1));
                bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                bias_ax.set_ylabel('Model Error, mm',{'fontsize': plt_param.LABELfs});    
                bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                ###################################################################
                #insert data statistics as title
                BT = bias_ax.set_title('RMSE: N/A',{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                BT.set_position(plt_param.Bloc_left);
                ###################################################################
                pwv_fig.subplots_adjust(left = 0.125,right = 0.935,bottom = 0.155, top = 0.925)
                ###################################################################
                bias_ax.annotate(Date_start_label,xy = plt_param.metDATEpad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                bias_ax.annotate(Date_end_label,xy = plt_param.metDATE_end_pad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                #######################################################################
                # Add some extra info
                bias_ax.annotate(MDLgridinfo,xy = plt_param.GRIDpad_left_dual,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'left',\
                                fontsize = plt_param.TITLEfs)
                bias_ax.annotate(ATTRstr,xy = plt_param.ATTRpad_ctr_dual,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'center',\
                                fontsize = plt_param.TITLEfs)
                #######################################################################
                # add a logo
                imax = pwv_fig.add_axes(plt_param.metLOGOaxcoords)
                imax.set_axis_off()
                imax.imshow(logo, aspect="equal")
                #######################################################################
                #Save the Developed Figure
                FILout = os.path.join(IMGpth_model,FILpfx+MDLname+'_PWV_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                pwv_fig.savefig(FILout,dpi=300); plt.close('all');
                ###################################################################
    
                """ ###############################################################
                #Boundary Layer Plots (instruct to share x-axis)
                ############################################################### """
                pblh_fig, (pblh_ax,bias_ax) = plt.subplots(2,1,sharex = True, figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                         facecolor = 'w',edgecolor ='w',clear='True')
                ###################################################################
                # Plot the Atmospheric Boundary Layer Time Histories
                ##############################################
                #Observations
                pblh_ax.plot(MDLtime,MEASpblh_mu,'s',color='gray',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Obs ('+AVGlabel+' Avg)')
                pblh_ax.plot(MEAStime,MEASpblh, color='black',linewidth=plt_param.PLTlwidth,label='Obs')
                pblh_ax.plot(MEAStime,MEASpblh,'.',color='black',markersize=1,linestyle='None',label='Obs')
                pblh_ax.plot(MDLtime,MDLtime*np.nan,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname+MDLstr)
                pblh_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                handles, labels = pblh_ax.get_legend_handles_labels()
                order = [1,0,3]
                pblh_ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], \
                               loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                               ncol = 3,frameon = False,columnspacing = plt_param.LEGcol_space);    
                ###################################################################
                pblh_ax.set_ylim(PBLHmin,PBLHmax); pblh_ax.set_yticks(PBLHtick);
                pblh_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                pblh_ax.set_ylabel('Boundary Layer,\nm AGL',{'fontsize': plt_param.LABELfs});
                pblh_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                ###################################################################
                #insert station info as title
                T = pblh_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                T.set_position(plt_param.Tloc);
                
                """ ###############################################################
                # Plot the Model Boundary Layer Biases
                ##############################################################  """
                # initate 0/0 line
                bias_ax.plot([Xmin,Xmax],[0,0],color='black',linewidth=plt_param.PLTlwidth)
                ###################################################################
                bias_ax.plot(MDLtime,MDLtime*np.nan,'^-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Mdl - Obs')
                bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                bias_ax.legend(loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                             ncol = 1,frameon = False);  
                ###################################################################
                bias_ax.set_xticks(Xtick); bias_ax.set_xticklabels(Xtick_lab,fontsize=plt_param.TICKfs);
                bias_ax.set_xlim(Xmin,Xmax)
                ###################################################################
                bias_ax.set_xlabel('Hour, UTC',{'fontsize': plt_param.LABELfs});
                bias_ax.xaxis.set_label_coords(plt_param.Xlab_pad_dual[0],plt_param.Xlab_pad_dual[1]);
                ###################################################################
                # Denote that Model Data was Not Availabel
                bias_ax.text(np.nanmean([Xmin,Xmax]),0.25,MDLname+' Data Not Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
                ###################################################################
                #Define Bounds of Wind Speed - Model Bias to Determine Y-Bounds
                bias_ax.set_ylim(-3,3); bias_ax.set_yticks(np.arange(-3,3.1,1));
                bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                #insert data statistics as title
                BT = bias_ax.set_title('RMSE: N/A',{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                BT.set_position(plt_param.Bloc_left);
                ###################################################################
                bias_ax.set_ylabel('Model Error, m',{'fontsize': plt_param.LABELfs});    
                bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                ###################################################################
                pblh_fig.subplots_adjust(left = 0.125,right = 0.935,bottom = 0.155, top = 0.925)
                ###################################################################
                bias_ax.annotate(Date_start_label,xy = plt_param.metDATEpad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                bias_ax.annotate(Date_end_label,xy = plt_param.metDATE_end_pad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                #######################################################################
                # Add some extra info
                bias_ax.annotate(MDLgridinfo,xy = plt_param.GRIDpad_left_dual,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'left',\
                                fontsize = plt_param.TITLEfs)
                bias_ax.annotate(ATTRstr,xy = plt_param.ATTRpad_ctr_dual,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'center',\
                                fontsize = plt_param.TITLEfs)
                #######################################################################
                # add a logo
                imax = pblh_fig.add_axes(plt_param.metLOGOaxcoords)
                imax.set_axis_off()
                imax.imshow(logo, aspect="equal")
                #######################################################################
                #Save the Developed Figure
                FILout = os.path.join(IMGpth_model,FILpfx+MDLname+'_PBLH_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                pblh_fig.savefig(FILout,dpi=300); plt.close('all');
                ###################################################################
    
                """ ###############################################################
                #Liquid Water Path Plots (instruct to share x-axis)
                ############################################################### """
                lwp_fig, (lwp_ax,bias_ax) = plt.subplots(2,1,sharex = True, figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                         facecolor = 'w',edgecolor ='w',clear='True')
                ###################################################################
                # Plot the Liquid Water Path Time Histories
                ##############################################
                #Observations
                lwp_ax.plot(MDLtime,MEASlwp_mu,'s',color='gray',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Obs ('+AVGlabel+' Avg)')
                lwp_ax.plot(MEAStime,MEASlwp, color='black',linewidth=plt_param.PLTlwidth,label='Obs')
                lwp_ax.plot(MEAStime,MEASlwp,'.',color='black',markersize=1,linestyle='None',label='Obs')
                lwp_ax.plot(MDLtime,MDLtime*np.nan,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname+MDLstr)
                lwp_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                handles, labels = lwp_ax.get_legend_handles_labels()
                order = [1,0,3]
                lwp_ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], \
                               loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                               ncol = 3,frameon = False,columnspacing = plt_param.LEGcol_space);    
                ###################################################################
                lwp_ax.set_ylim(LWPmin,LWPmax); lwp_ax.set_yticks(LWPtick);
                lwp_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                lwp_ax.set_ylabel('Liquid Water Path,\nm m-2',{'fontsize': plt_param.LABELfs});
                lwp_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                ###################################################################
                #insert station info as title
                T = lwp_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                T.set_position(plt_param.Tloc);
                
                """ ###############################################################
                # Plot the Model Liquid Water Path Biases
                ##############################################################  """
                # initate 0/0 line
                bias_ax.plot([Xmin,Xmax],[0,0],color='black',linewidth=plt_param.PLTlwidth)
                ###################################################################
                bias_ax.plot(MDLtime,MDLtime*np.nan,'^-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Mdl - Obs')
                bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                bias_ax.legend(loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                             ncol = 1,frameon = False);  
                ###################################################################
                bias_ax.set_xticks(Xtick); bias_ax.set_xticklabels(Xtick_lab,fontsize=plt_param.TICKfs);
                bias_ax.set_xlim(Xmin,Xmax)
                ###################################################################
                bias_ax.set_xlabel('Hour, UTC',{'fontsize': plt_param.LABELfs});
                bias_ax.xaxis.set_label_coords(plt_param.Xlab_pad_dual[0],plt_param.Xlab_pad_dual[1]);
                ###################################################################
                # Denote that Model Data was Not Availabel
                bias_ax.text(np.nanmean([Xmin,Xmax]),0.25,MDLname+' Data Not Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
                ###################################################################
                #Define Bounds of Wind Speed - Model Bias to Determine Y-Bounds
                bias_ax.set_ylim(-3,3); bias_ax.set_yticks(np.arange(-3,3.1,1));
                bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                #insert data statistics as title
                BT = bias_ax.set_title('RMSE: N/A',{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                BT.set_position(plt_param.Bloc_left);
                ###################################################################
                bias_ax.set_ylabel('Model Error, g m-2',{'fontsize': plt_param.LABELfs});    
                bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                ###################################################################
                lwp_fig.subplots_adjust(left = 0.125,right = 0.935,bottom = 0.155, top = 0.925)
                ###################################################################
                bias_ax.annotate(Date_start_label,xy = plt_param.metDATEpad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                bias_ax.annotate(Date_end_label,xy = plt_param.metDATE_end_pad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                #######################################################################
                # Add some extra info
                bias_ax.annotate(MDLgridinfo,xy = plt_param.GRIDpad_left_dual,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'left',\
                                fontsize = plt_param.TITLEfs)
                bias_ax.annotate(ATTRstr,xy = plt_param.ATTRpad_ctr_dual,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'center',\
                                fontsize = plt_param.TITLEfs)
                #######################################################################
                # add a logo
                imax = lwp_fig.add_axes(plt_param.metLOGOaxcoords)
                imax.set_axis_off()
                imax.imshow(logo, aspect="equal")
                #######################################################################
                #Save the Developed Figure
                FILout = os.path.join(IMGpth_model,FILpfx+MDLname+'_LWP'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                lwp_fig.savefig(FILout,dpi=300); plt.close('all');
                ###################################################################
