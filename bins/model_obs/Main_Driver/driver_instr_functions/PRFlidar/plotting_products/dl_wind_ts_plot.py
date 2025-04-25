#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This Function Plots for the Case where both Measurement and Model Data are
Available
    
###############################################################################
Created on Tue, Feb 27, 2024
@author: dmurray after jduncan
###############################################################################
The Purpose of this Script is to Develop and Output Wind-Related Plots Based on
Measurements from the Radar Wind Profiler Systems Installed at 
Various Sites
###############################################################################
Outputted Plots Follow:
instrument_whatisplotted_Sid.date.png
**Imperative that the "." is the first period -- this is used to strip the **
**date from the file so that duplicate images are not run **
###############################################################################
Most Plots Generated in This Script Will Consist of Two Panels 
- in the first plot, the obs will be plotted along with the model output
- in the second plot, the model bias will be shown and statistics will be annotated

- Model Data Still Needs to be Incorporated
###############################################################################
    #For Statisical Wind Direction Computations -- Only Consider Wind Speeds 
    #Abovce Some Threshold Value

"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as img
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection
import os,time
from datetime import datetime,timedelta
#####################################
from coord_distance import dist2coord
from atmo_meas_to_hour_end  import atmo_meas2hr_end
from atmo_meas_to_hour_mid  import atmo_meas2hr_mid
from wd_meas_to_hour_end import wd_meas2hr_end
from wd_meas_to_hour_mid import wd_meas2hr_mid
from time_tuple_cat import time_cat
from atmo_time_correction import atmo_time_correct
from atmo_spatial_tuple_cat import atmo_spatial_cat
from ws_to_power import speed2power_offshore
from ws_to_power import ramp_up_down_segs
#####################################
from atmo_bias import atmo_bias
from wd_bias_flip import wd_bias
#####################################
from wd_ylim_plot import wd_ylim
from wd_ylim_bias import wd_blim
from wd_disc import wd_plt
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
AVGper = getattr(plt_param,'DL_avgmins',30)
AVGlabel = avgmins_to_label(AVGper)
###############################################################################

def dl_WTSplot(MEASin,MDLin,OUTinfo):
    DATE = OUTinfo[3];
    basedate=datetime.strptime(DATE,'%Y%m%d')
    ###########################################################################
    IMGpth = os.path.join(OUTinfo[0], DATE);
    if not os.path.exists(IMGpth):
        os.makedirs(IMGpth)
    IMGpth_ext = os.path.join(OUTinfo[0],'ext', DATE);
    if not os.path.exists(IMGpth_ext):
        os.makedirs(IMGpth_ext)
    IMGpth_model = os.path.join(OUTinfo[0],'no_model', DATE);
    if not os.path.exists(IMGpth_model):
        os.makedirs(IMGpth_model)
    ###########################################################################
    MDLname = OUTinfo[1];
    Sid = OUTinfo[2];
    logo=img.imread(plt_param.NOAAlogo)
    FILpfx='dlwinds_'
    ATTRstr = "Credit: University of Colorado Boulder/ATOC"
    
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
    psl_WTSplot(MEASin,MDLin,FILout)
    #############################################
    OUTinfo :: FILout = (IMGout,MDLname,Sid,doi)
    #############################################"""
    ###########################################################################
    """ Unpack Modeled Wind Data
    ###########################
    MDLin :: 
            HRRRin = (ESRLhrrr_ini,ESRLhrrr_loc,ESRLhrrr_xy,ESRLhrrr_wind)
            RAPin = (ESRLrap_ini,ESRLrap_loc,ESRLrap_xy,ESRLrap_wind)
            #-> each embedded tuple has length of initialization files
            ###################################################################
            MDLloc = (PSLsite,PSLlon,PSLlat,Mlon,Mlat,Melev,Mdist)
            MDLxy = (Mfcst_tme,Mmeas_hgt)
            MDLwind = depends on interpolation choice (see below)
                    ESRLw_o = (WSo,WDo,WSsfc_o,WDsfc_o,Uo,Usfc_o,Vo,Vsfc_o,Wo)
                    ESRLw_int = (WSint,WDint,WDsfc_o,WDsfc_int,Uint,Usfc_int,Vint,Vsfc_int,Wint)   
                    ###########################################################################
                    ESRLw = (ESRLw_o,ESRLw_int)
    ####################################################################### """
    # Model Initialization Times #        
    MDLini = MDLin[0]; #outputs model initialization time for each ini hour
    MDLxy = MDLin[2]; #outputs forecast time and measurement height
    MDLwind = MDLin[3] #outputs wind ini hour wind data (expand in loop)
    ###########################################################################
    MDLloc = MDLin[1][0]; #just store first initialization hour
    MDLloc_names = MDLloc[0];
    MDLlon_sites = MDLloc[1]; MDLlat_sites = MDLloc[2]; # desired grid point
    ######################################################################
    MDLlon_grids = MDLloc[3]; MDLlat_grids = MDLloc[4];# model grid point referenced
    MDLelev_grids = MDLloc[5]
    MDLdist2sites = MDLloc[6]

    ###########################################################################
    """ Define Measurement Location and Model Grid Point Structure """
    MEASlon = getattr(station_param,Sid+'lon')
    MEASlat = getattr(station_param,Sid+'lat')
    MEASalt = getattr(station_param,Sid+'alt')
    MEASid = getattr(station_param,Sid+'name')
    MEASname = getattr(station_param,Sid+'descr')
    HUBhgt = getattr(station_param,Sid+'hubhgt')
    ###########################################################################
    # Get alternate grid point if available
    MEASmdlid = getattr(station_param,Sid+'mdlID',"") # mdl ID for multiple points
    ALTsite = getattr(station_param,Sid+'altSite',None)
    ALTmdlid = getattr(station_param,Sid+'altmdlID',"") # mdl ID for alternate point
    ALTmdls = getattr(station_param,Sid+'altmdls',[]) # mdl ID for alternate point
    #MEAStitle = MEASname+' ('+str(MEASalt)+' m)'
    MEAStitle = MEASname+'\n'+format_lat_lon_alt(MEASlat,MEASlon,MEASalt)
    ###########################################################################
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
    """ Unpack Measured Wind Information 
    ####################################
    MEASin :: DLw_plt = (DLloc,DLtme,DLw_in)   
    ---> DLloc = (DLlon,DLlat,DLalt)
    ---> DLtme = (DLxy_minus,DLxy,DLxy_plus)
    ---> DLw_in = (DLwind_minus,DLwind,DLwind_plus)
    ###########################################################################
    DLwind = (WS,WD,Uvel,Vvel,Wvel,SNR)
    ########################################################################"""
    # Declare Measurement Heights
    MEAShgt = MEASin[1][1][1];
    
    # Define Relevant Indices for Reference in Concatenation #
    WSind_cat = 0; WDind_cat = 1;
    Uind_cat = 2; Vind_cat = 3;
    ##########################################
    # SNRind_cat = 5; #non-range-corrected SNR
    SNRind_cat = 6; #range-corected SNR
    ###########################################################################
    # For Spatial Concatenation -- Define the Time Axis within the Input Tuple
    Tax = 1;
        
    #-> Combine Relevant Measurement Tuples 
    MEAStime = time_cat(MEASin[1][0],MEASin[1][1],MEASin[1][2],MEASin[1][3])
    ###########################################################################
    MEASw_in = MEASin[2]; 
    #####################
    MEASws = atmo_spatial_cat(MEASw_in[0],MEASw_in[1],MEASw_in[2],MEASw_in[3],WSind_cat,Tax)
    MEASwd = atmo_spatial_cat(MEASw_in[0],MEASw_in[1],MEASw_in[2],MEASw_in[3],WDind_cat,Tax)
    ################################################################
    
    MEAShubidx = np.argmin(np.abs(MEAShgt-HUBhgt));
    MEAShubhgt_value = int(MEAShgt[MEAShubidx])

    ###########################################################################
    #  Extract out the hub data from the low level winds
    MEASws = MEASws[MEAShubidx,:]
    MEASwd = MEASwd[MEAShubidx,:]

    ###################################################################
    # Perform Correction to Create Equally Spaced Data Arrays #
    ###################################################################
    MEAStime_new, MEASws = atmo_time_correct(MEAStime,MEASws);
    _, MEASwd = atmo_time_correct(MEAStime,MEASwd);
    MEAStime = MEAStime_new

    MEASpow = speed2power_offshore(MEASws)

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
            # Determine Whether Model Data Exists for Said Initialization Time
            if INIpos[ini] in MDLini:
                """ Model Initialization Time Defined """
                ###################################################################
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
                MDLhgt = MDLxy[INIref][1]; #model height information
                MDLhgt_ref = np.abs(MDLhgt-HUBhgt).argmin();
                MDLhgt_ref_value = int(MDLhgt[MDLhgt_ref])
                ####################################################
                MDLw_ini = MDLwind[INIref][INTPoi]
                ###################################################################
                # gridded wind information (WSo,WDo,WSsfc_o,WDsfc_o,Uo,Usfc_o,Vo,Vsfc_o,Wo)
                MDLws = MDLw_ini[0]; MDLwd = MDLw_ini[1];
                # data is currently shaped as: (forecast_hour,sites,hgt)
                ###################################################################
                # extract from relevant grid point and level
                ###################################################################
                MDLws = MDLws[:,Gref,MDLhgt_ref]
                MDLwd = MDLwd[:,Gref,MDLhgt_ref]
                MDLpow = speed2power_offshore(MDLws)
                MDLws_up_segs, MDLws_down_segs = ramp_up_down_segs(MDLtime,MDLws,threshold=plt_param.WSthreshold)
                MDLpow_up_segs, MDLpow_down_segs = ramp_up_down_segs(MDLtime,MDLpow,threshold=plt_param.POWthreshold)
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
                    MDLws_ext = MDLws; MDLwd_ext = MDLwd; MDLpow_ext = MDLpow;
                    MDLtime_ext = MDLtime;
                    MDLws_ext_up_segs, MDLws_ext_down_segs = ramp_up_down_segs(MDLtime,MDLws_ext,threshold=plt_param.WSthreshold)
                    MDLpow_ext_up_segs, MDLpow_ext_down_segs = ramp_up_down_segs(MDLtime,MDLpow_ext,threshold=plt_param.POWthreshold)
                    ###############################################################
                    # Perform Model Averaging for Plotting and Statistical Analsyes
                    MEASws_ext_mu = atmo_meas2hr_mid(MEASws,None,MEAStime,MDLtime_ext,AVGmins=AVGper,func='median');
                    MEASpow_ext_mu = atmo_meas2hr_mid(MEASpow,None,MEAStime,MDLtime_ext,AVGmins=AVGper,func='median');
                    MEASwd_ext_mu = wd_meas2hr_mid(MEASwd,None,MEAStime,MDLtime_ext,AVGmins=AVGper,func='median');
                    MEASws_ext_up_segs, MEASws_ext_down_segs = ramp_up_down_segs(MDLtime,MEASws_ext_mu,threshold=plt_param.WSthreshold)
                    MEASpow_ext_up_segs, MEASpow_ext_down_segs = ramp_up_down_segs(MDLtime,MEASpow_ext_mu,threshold=plt_param.POWthreshold)
                    ###############################################################
                    WSmin_ext,WSmax_ext,WStick_ext = atmo_bound_lim((MEASws[ext_meas_oi],\
                                                   MDLws),plt_param.WSrnd_bse,plt_param.WSmin,plt_param.WSmax,\
                                                   plt_param.WSrnd_bse)
                    #POWmin_ext,POWmax_ext,POWtick_ext = atmo_bound_lim((MEASpow[ext_meas_oi],\
                    #                               MDLpow),plt_param.POWrnd_bse,plt_param.POWmin,plt_param.POWmax,\
                    #                               plt_param.POWrnd_bse)
                    POWmin_ext=plt_param.POWmin;POWmax_ext=plt_param.POWmax;POWtick_ext=np.arange(0,plt_param.POWmax,plt_param.POWrnd_bse);
                    #WDmin_ext,WDmax_ext,WDtick_ext = wd_ylim((MEASwd[ext_meas_oi],MDLwd_ext));
                    WDmin_ext=0;WDmax_ext=360;WDtick_ext=np.arange(0,361,60);
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
                    
                    """ Define Bounds for the Model-Specific Standard Forecast Period """
                    meas_oi = np.logical_and(MEAStime >= MDLtime[0], MEAStime <= MDLtime[FXlen-1]) 
                    ###############################################################
                    # Adjust Model Variables #
                    MDLws = MDLws[0:FXlen]; MDLwd = MDLwd[0:FXlen]; MDLpow = MDLpow[0:FXlen];
                    MDLtime = MDLtime[0:FXlen];
                    MDLws_up_segs, MDLws_down_segs = ramp_up_down_segs(MDLtime,MDLws,threshold=plt_param.WSthreshold)
                    MDLpow_up_segs, MDLpow_down_segs = ramp_up_down_segs(MDLtime,MDLpow,threshold=plt_param.POWthreshold)
                    ###############################################################
                    ###############################################################
                    # Perform Model Averaging for Plotting and Statistical Analsyes
                    MEASws_mu = atmo_meas2hr_mid(MEASws,None,MEAStime,MDLtime,AVGmins=AVGper,func='median');
                    MEASwd_mu = wd_meas2hr_mid(MEASwd,None,MEAStime,MDLtime,AVGmins=AVGper,func='median');
                    MEASpow_mu = atmo_meas2hr_mid(MEASpow,None,MEAStime,MDLtime,AVGmins=AVGper,func='median');
                    MEASws_up_segs, MEASws_down_segs = ramp_up_down_segs(MDLtime,MEASws_mu,threshold=plt_param.WSthreshold)
                    MEASpow_up_segs, MEASpow_down_segs = ramp_up_down_segs(MDLtime,MEASpow_mu,threshold=plt_param.POWthreshold)
                    ###############################################################
                    WSmin,WSmax,WStick = atmo_bound_lim((MEASws[meas_oi],\
                                                   MDLws),plt_param.WSrnd_bse,plt_param.WSmin,plt_param.WSmax,\
                                                   plt_param.WSrnd_bse)
                    #POWmin,POWmax,POWtick = atmo_bound_lim((MEASpow[meas_oi],\
                    #                               MDLpow),plt_param.POWrnd_bse,plt_param.POWmin,plt_param.POWmax,\
                    #                               plt_param.POWrnd_bse)
                    POWmin=plt_param.POWmin;POWmax=plt_param.POWmax;POWtick=np.arange(0,plt_param.POWmax,plt_param.POWrnd_bse);
                    #WDmin,WDmax,WDtick = wd_ylim((MEASwd[meas_oi],MDLwd));
                    WDmin=0;WDmax=360;WDtick=np.arange(0,361,60);
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
                    # no model adjustment needed 
                    # or less due to a download issue (i.e. MDLws = MDLws[0:FXlen])
                    ###############################################################
                    # Perform Model Averaging for Plotting and Statistical Analsyes
                    MEASws_mu = atmo_meas2hr_mid(MEASws,None,MEAStime,MDLtime,AVGmins=AVGper,func='median');
                    MEASwd_mu = wd_meas2hr_mid(MEASwd,None,MEAStime,MDLtime,AVGmins=AVGper,func='median');
                    MEASpow_mu = atmo_meas2hr_mid(MEASpow,None,MEAStime,MDLtime,AVGmins=AVGper,func='median');
                    MEASws_up_segs, MEASws_down_segs = ramp_up_down_segs(MDLtime,MEASws_mu,threshold=plt_param.WSthreshold)
                    MEASpow_up_segs, MEASpow_down_segs = ramp_up_down_segs(MDLtime,MEASpow_mu,threshold=plt_param.POWthreshold)
                    # desire mean values at each desired forecast time (push forward)
                    ###############################################################
                    WSmin,WSmax,WStick = atmo_bound_lim((MEASws[meas_oi],\
                                                   MDLws),plt_param.WSrnd_bse,plt_param.WSmin,plt_param.WSmax,\
                                                   plt_param.WSrnd_bse)
                    #POWmin,POWmax,POWtick = atmo_bound_lim((MEASpow[meas_oi],\
                    #                               MDLpow),plt_param.POWrnd_bse,plt_param.POWmin,plt_param.POWmax,\
                    #                               plt_param.POWrnd_bse)
                    POWmin=plt_param.POWmin;POWmax=plt_param.POWmax;POWtick=np.arange(0,plt_param.POWmax,plt_param.POWrnd_bse);
                    #WDmin,WDmax,WDtick = wd_ylim((MEASwd[meas_oi],MDLwd));
                    WDmin=0;WDmax=360;WDtick=np.arange(0,361,60);
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
                # would not enter loop if initialization time was not available
                """ ###############################################################
                Wind Speed Plots (instruct to share x-axis)
                ############################################################### """
                ws_fig, (ws_ax,bias_ax) = plt.subplots(2,1,sharex = True, figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                            facecolor = 'w',edgecolor ='w',clear='True')
                ###################################################################
                # Plot the Wind Speed Time Histories
                ####################################
                #Observations -- Only Plotting the QC'd Values
                ws_ax.plot(MDLtime,MEASws_mu,'s-',color='gray',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Obs ('+AVGlabel+' Avg)')
                ws_ax.plot(MEAStime,MEASws,'.',color='black',markersize=1,linewidth=0,label='Obs')
                ws_ax.plot(MDLtime,MDLws,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname+' ('+str(MDLhgt_ref_value)+' m)')
                ws_ax.add_collection(LineCollection(MEASws_up_segs,colors=(plt_param.WSramp_up_color),linewidths=(plt_param.PLTlwidth),label='Up ramps'))
                ws_ax.add_collection(LineCollection(MEASws_down_segs,colors=(plt_param.WSramp_down_color),linewidths=(plt_param.PLTlwidth),label='Down ramps'))
                ws_ax.add_collection(LineCollection(MDLws_up_segs,colors=(plt_param.WSramp_up_color),linewidths=(plt_param.PLTlwidth),label=' '))
                ws_ax.add_collection(LineCollection(MDLws_down_segs,colors=(plt_param.WSramp_down_color),linewidths=(plt_param.PLTlwidth),label=' '))
                ws_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                handles, labels = ws_ax.get_legend_handles_labels()
                # override line label 5 to be a dummy white line
                handles[5] = Line2D([0],[0],color="w")
                order = [1,5,0,3,2,4]
                ws_ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], \
                               loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                               ncol = 3,frameon = False,columnspacing = plt_param.LEGcol_space);    
                ###################################################################
                ws_ax.set_ylim(WSmin,WSmax); ws_ax.set_yticks(WStick);
                ws_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                ws_ax.set_ylabel('Wind Speed, m $\mathrm{s}^{-1}$',{'fontsize': plt_param.LABELfs});
                ws_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                ###################################################################
                #insert data date as title
                T = ws_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                T.set_position(plt_param.Tloc);
                
                """ ###########################################################
                Define and Plot the Model Wind Speed Biases
                ############################################################"""
                # initate 0/0 line
                bias_ax.plot([Xmin,Xmax],[0,0],color='black',linewidth=1.5)
                ###################################################################
                # Define the Bias (Observations - Measurements)
                WSbias = MDLws - MEASws_mu
                bias_ax.plot(MDLtime,WSbias,'^-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Mdl - Obs')
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
                #Define Bounds of Wind Speed - Model Error to Determine Y-Bounds
                WSb_min,WSb_max,WSb_tick = atmo_bound_lim((WSbias),plt_param.WSbias_bse,plt_param.WSbmin,\
                                                     plt_param.WSbmax,plt_param.WSbias_bse)
                bias_ax.set_ylim(WSb_min,WSb_max); bias_ax.set_yticks(WSb_tick);
                bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                bias_ax.set_ylabel('Model Error',{'fontsize': plt_param.LABELfs});    
                bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                ###################################################################
                #input basic statistics
                RMSEws,BIASws = atmo_bias(WSbias)
                #insert data statistics as title
                BT = bias_ax.set_title('RMSE: ' + str(RMSEws)+' m '+r'$\mathrm{s}^{-1}$ | Mean Bias: '+\
                                      str(BIASws)+' m '+r'$\mathrm{s}^{-1}$',\
                                      {'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5) 
                BT.set_position(plt_param.Bloc_left);
                ###################################################################
                #-> annotate measurement heights
                ws_ax.annotate(str(MEAShubhgt_value)+' m AGL',xy=(0.13,0.895), xycoords='figure fraction',color='black',\
                                 fontsize=plt_param.LEGfs,weight='bold')
                ###################################################################
                ws_fig.subplots_adjust(left = 0.125,right = 0.935,bottom = 0.155, top = 0.925)
                ###################################################################
                bias_ax.annotate(Date_start_label,xy = plt_param.metDATEpad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                                #fontsize = plt_param.IDfs,bbox = dict(boxstyle='square',fc='w'))
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
                imax = ws_fig.add_axes(plt_param.metLOGOaxcoords)
                imax.set_axis_off()
                imax.imshow(logo, aspect="equal")
                ###################################################################
                #Save the Developed Figure
                FILout = os.path.join(IMGpth,FILpfx+MDLname+'_WShub_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                ws_fig.savefig(FILout,dpi=300); plt.close('all');
                
                """ ###############################################################
                Wind Direction Plots (instruct to share x-axis)
                ############################################################### """
                wd_fig, (wd_ax,bias_ax) = plt.subplots(2,1,sharex = True, figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                         facecolor = 'w',edgecolor ='w',clear='True')
                ###################################################################
                # Plot the Wind Direction Time Histories
                ###################################################################
                wd_ax = wd_plt(MDLtime,MEASwd_mu,wd_ax,'s-','gray','Obs ('+AVGlabel+' Avg)',marker_size=plt_param.PLTmsize)
                wd_ax = wd_plt(MEAStime,MEASwd,wd_ax,'.','black','Obs',marker_size=1,line_width=0)
                wd_ax = wd_plt(MDLtime,MDLwd,wd_ax,'o-','red',MDLname+' ('+str(MDLhgt_ref_value)+' m)',marker_size=plt_param.PLTmsize,line_width=plt_param.PLTlwidth)
                ###################################################################
                wd_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                handles, labels = wd_ax.get_legend_handles_labels()
                order = [1,0,2]
                wd_ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], \
                               loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                               ncol = 3,frameon = False,columnspacing = plt_param.LEGcol_space);    
                ###################################################################
                wd_ax.set_ylim(WDmin,WDmax); wd_ax.set_yticks(WDtick);
                wd_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                wd_ax.set_ylabel('Wind Direction, deg ${}^{\circ}$',{'fontsize': plt_param.LABELfs});
                wd_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                ###################################################################
                #insert data date as title
                T = wd_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                T.set_position(plt_param.Tloc);
                
                """ ###############################################################
                # Plot the Model Wind Direction Biases
                ############################################################### """
                # initate 0/0 line
                bias_ax.plot([Xmin,Xmax],[0,0],color='black',linewidth=1.5)
                ###################################################################
                # Define the Bias (Observations - Measurements)
                WDbias,RMSEwd,BIASwd = wd_bias(MEASwd_mu,MDLwd)
                bias_ax.plot(MDLtime,WDbias,'^-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Mdl - Obs')
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
                #Define Bounds of Wind Speed - Model Bias to Determine Y-Bounds
                WDb_min,WDb_max,WDb_tick = wd_blim((WDbias));
                bias_ax.set_ylim(WDb_min,WDb_max); bias_ax.set_yticks(WDb_tick);
                bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                bias_ax.set_ylabel('Model Error',{'fontsize': plt_param.LABELfs});    
                bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                ###################################################################
                #insert data statistics as title
                BT = bias_ax.set_title('RMSE: ' + str(RMSEwd)+' deg '+r'$\mathrm{}^{\circ}$ | Mean Bias: '+\
                                      str(BIASwd)+' deg '+r'$\mathrm{}^{\circ}$',\
                                      {'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                BT.set_position(plt_param.Bloc_left);
                ###################################################################
                #-> annotate measurement heights
                wd_ax.annotate(str(MEAShubhgt_value)+' m AGL',xy=(0.13,0.895), xycoords='figure fraction',color='black',\
                                 fontsize=plt_param.LEGfs,weight='bold')
                ###################################################################
                wd_fig.subplots_adjust(left = 0.125,right = 0.935,bottom = 0.155, top = 0.925)
                ###################################################################
                bias_ax.annotate(Date_start_label,xy = plt_param.metDATEpad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                                #fontsize = plt_param.IDfs,bbox = dict(boxstyle='square',fc='w'))
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
                imax = wd_fig.add_axes(plt_param.metLOGOaxcoords)
                imax.set_axis_off()
                imax.imshow(logo, aspect="equal")
                ###################################################################
                #Save the Developed Figure
                FILout = os.path.join(IMGpth,FILpfx+MDLname+'_WDhub_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                wd_fig.savefig(FILout,dpi=300); plt.close('all')
                ###################################################################
                
                """ ###############################################################
                Hub Power Profile (instruct to share x-axis)
                ############################################################### """
                pow_fig, (pow_ax,bias_ax) = plt.subplots(2,1,sharex = True, figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                         facecolor = 'w',edgecolor ='w',clear='True')
                ####################################################################
                # Plot the Power Time Histories
                #####################################
                #Observations -- Only Plotting the QC'd Values
                pow_ax.plot(MDLtime,MEASpow_mu,'s-',color='gray',markersize=plt_param.PLTmsize,label='Obs ('+AVGlabel+' Avg)')
                pow_ax.plot(MEAStime,MEASpow,'.',color='black',markersize=1,linewidth=0,label='Obs')
                pow_ax.plot(MDLtime,MDLpow,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname+' ('+str(MDLhgt_ref_value)+' m)')
                pow_ax.add_collection(LineCollection(MEASpow_up_segs,colors=(plt_param.POWramp_up_color),linewidths=(plt_param.PLTlwidth),label='Up ramps'))
                pow_ax.add_collection(LineCollection(MEASpow_down_segs,colors=(plt_param.POWramp_down_color),linewidths=(plt_param.PLTlwidth),label='Down ramps'))
                pow_ax.add_collection(LineCollection(MDLpow_up_segs,colors=(plt_param.POWramp_up_color),linewidths=(plt_param.PLTlwidth),label=' '))
                pow_ax.add_collection(LineCollection(MDLpow_down_segs,colors=(plt_param.POWramp_down_color),linewidths=(plt_param.PLTlwidth),label=' '))
                pow_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                handles, labels = pow_ax.get_legend_handles_labels()
                # override line label 5 to be a dummy white line
                handles[5] = Line2D([0],[0],color="w")
                order = [1,5,0,3,2,4]
                pow_ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], \
                               loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                               ncol = 3,frameon = False,columnspacing = plt_param.LEGcol_space);    
                ###################################################################
                pow_ax.set_ylim(POWmin,POWmax); pow_ax.set_yticks(POWtick);
                pow_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                pow_ax.set_ylabel('Wind Capacity Factor',{'fontsize': plt_param.LABELfs});
                pow_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                ###################################################################
                #insert station info as title
                T = pow_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                T.set_position(plt_param.Tloc);
                             
                """ ###############################################################
                # Plot the Model Power Biases
                ############################################################### """
                # initate 0/0 line
                bias_ax.plot([Xmin,Xmax],[0,0],color='black',linewidth=1.5)
                ###################################################################
                # Define the Bias (Observations - Measurements)
                POWbias = MDLpow - MEASpow_mu
                bias_ax.plot(MDLtime,POWbias,'^-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Mdl - Obs')
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
                #Define Bounds of Power - Model Bias to Determine Y-Bounds
                #POWb_min,POWb_max,POWb_tick = atmo_bound_lim((POWbias),plt_param.POWbias_bse,plt_param.POWbmin,\
                #                                     plt_param.POWbmax,plt_param.POWbias_bse)
                POWb_min=-1.1;POWb_max=1.1;POWb_tick=np.arange(-1,1.1,.5);
                bias_ax.set_ylim(POWb_min,POWb_max); bias_ax.set_yticks(POWb_tick);
                bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                bias_ax.set_ylabel('Model Error',{'fontsize': plt_param.LABELfs});    
                bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                ###################################################################
                #input basic statistics
                RMSEpow,BIASpow = atmo_bias(POWbias)
                #insert data statistics as title
                BT = bias_ax.set_title('RMSE: ' + str(RMSEpow)+' | Mean Bias: '+\
                                      str(BIASpow)+' ',\
                                      {'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5) 
                BT.set_position(plt_param.Bloc_left);
                ###################################################################
                #-> annotate measurement heights
                pow_ax.annotate(str(MEAShubhgt_value)+' m AGL',xy=(0.13,0.895), xycoords='figure fraction',color='black',\
                                 fontsize=plt_param.LEGfs,weight='bold')
                #######################################################################
                pow_fig.subplots_adjust(left = 0.125,right = 0.935,bottom = 0.155, top = 0.925)
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
                imax = pow_fig.add_axes(plt_param.metLOGOaxcoords)
                imax.set_axis_off()
                imax.imshow(logo, aspect="equal")
                #######################################################################
                #Save the Developed Figure
                FILout = os.path.join(IMGpth,FILpfx+MDLname+'_POW_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                pow_fig.savefig(FILout,dpi=300); plt.close('all');
                """ Plot Extended Forecast Period if Applicable """
                ###################################################################
                if 'MDLtime_ext' in locals():
                    """ Extended Wind Speed Forecast """
                    ws_fig, (ws_ax,bias_ax) = plt.subplots(2,1,sharex = True, figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                            facecolor = 'w',edgecolor ='w',clear='True')
                    ###############################################################
                    # Plot the Wind Speed Time Histories
                    ####################################
                    #Observations -- Only Plotting the QC'd Values
                    ws_ax.plot(MDLtime_ext,MEASws_ext_mu,'s-',color='gray',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Obs ('+AVGlabel+' Avg)')
                    ws_ax.plot(MEAStime,MEASws,'.',color='black',markersize=1,linewidth=0,label='Obs')
                    ws_ax.plot(MDLtime_ext,MDLws_ext,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname+' ('+str(MDLhgt_ref_value)+' m)')
                    ws_ax.add_collection(LineCollection(MEASws_ext_up_segs,colors=(plt_param.WSramp_up_color),linewidths=(plt_param.PLTlwidth),label='Up ramps'))
                    ws_ax.add_collection(LineCollection(MEASws_ext_down_segs,colors=(plt_param.WSramp_down_color),linewidths=(plt_param.PLTlwidth),label='Down ramps'))
                    ws_ax.add_collection(LineCollection(MDLws_ext_up_segs,colors=(plt_param.WSramp_up_color),linewidths=(plt_param.PLTlwidth),label=' '))
                    ws_ax.add_collection(LineCollection(MDLws_ext_down_segs,colors=(plt_param.WSramp_down_color),linewidths=(plt_param.PLTlwidth),label=' '))
                    ws_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                    ###################################################################
                    handles, labels = ws_ax.get_legend_handles_labels()
                    # override line label 5 to be a dummy white line
                    handles[5] = Line2D([0],[0],color="w")
                    order = [1,5,0,3,2,4]
                    ws_ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], \
                                   loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                                   ncol = 3,frameon = False,columnspacing = plt_param.LEGcol_space);    
                    ###############################################################
                    ws_ax.set_ylim(WSmin_ext,WSmax_ext); ws_ax.set_yticks(WStick_ext);
                    ws_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                    ###############################################################
                    ws_ax.set_ylabel('Wind Speed, m $\mathrm{s}^{-1}$',{'fontsize': plt_param.LABELfs});
                    ws_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                    ###############################################################
                    #insert data date as title
                    T = ws_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                    T.set_position(plt_param.Tloc);
                    
                    """ ###########################################################
                    Define and Plot the Model Wind Speed Biases
                    ########################################################### """
                    # initate 0/0 line
                    bias_ax.plot([Xmin_ext,Xmax_ext],[0,0],color='black',linewidth=1.5)
                    ###############################################################
                    # Define the Bias (Observations - Measurements)
                    WSbias = MDLws_ext - MEASws_ext_mu
                    bias_ax.plot(MDLtime_ext,WSbias,'^-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Mdl - Obs')
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
                    #Define Bounds of Wind Speed - Model Error to Determine Y-Bounds
                    WSb_min,WSb_max,WSb_tick = atmo_bound_lim((WSbias),plt_param.WSbias_bse,plt_param.WSbmin,\
                                                     plt_param.WSbmax,plt_param.WSbias_bse)
                    bias_ax.set_ylim(WSb_min,WSb_max); bias_ax.set_yticks(WSb_tick);
                    bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                    ###############################################################
                    bias_ax.set_ylabel('Model Error',{'fontsize': plt_param.LABELfs});    
                    bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                    ###############################################################
                    #input basic statistics
                    RMSEws,BIASws = atmo_bias(WSbias)
                    #insert data statistics as title
                    BT = bias_ax.set_title('RMSE: ' + str(RMSEws)+' m '+r'$\mathrm{s}^{-1}$ | Mean Bias: '+\
                                          str(BIASws)+' m '+r'$\mathrm{s}^{-1}$',\
                                          {'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5) 
                    BT.set_position(plt_param.Bloc_left);
                    ###############################################################
                    #-> annotate measurement heights
                    ws_ax.annotate(str(MEAShubhgt_value)+' m AGL',xy=(0.13,0.895), xycoords='figure fraction',color='black',\
                                     fontsize=plt_param.LEGfs,weight='bold')
                    ###############################################################
                    ws_fig.subplots_adjust(left = 0.125,right = 0.935,bottom = 0.155, top = 0.925)
                    ###############################################################
                    bias_ax.annotate(Date_start_label_ext,xy = plt_param.metDATEpad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                                #fontsize = plt_param.IDfs,bbox = dict(boxstyle='square',fc='w'))
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
                    imax = ws_fig.add_axes(plt_param.metLOGOaxcoords)
                    imax.set_axis_off()
                    imax.imshow(logo, aspect="equal")
                    ###################################################################
                    #Save the Developed Figure
                    FILout = os.path.join(IMGpth_ext,FILpfx+MDLname+'_WShub_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                    ws_fig.savefig(FILout,dpi=300); plt.close('all');
                    
                    """ ###########################################################
                    Wind Direction Plots (instruct to share x-axis)
                    ########################################################### """
                    wd_fig, (wd_ax,bias_ax) = plt.subplots(2,1,sharex = True, figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                             facecolor = 'w',edgecolor ='w',clear='True')
                    ###############################################################
                    # Plot the Wind Direction Time Histories
                    ###############################################################
                    wd_ax = wd_plt(MDLtime_ext,MEASwd_ext_mu,wd_ax,'s-','gray','Obs ('+AVGlabel+' Avg)',marker_size=plt_param.PLTmsize)
                    wd_ax = wd_plt(MEAStime,MEASwd,wd_ax,'.','black','Obs',marker_size=1,line_width=0)
                    wd_ax = wd_plt(MDLtime_ext,MDLwd_ext,wd_ax,'o-','red',MDLname+' ('+str(MDLhgt_ref_value)+' m)',marker_size=plt_param.PLTmsize,line_width=plt_param.PLTlwidth)
                    ###############################################################
                    wd_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                    ###############################################################
                    handles, labels = wd_ax.get_legend_handles_labels()
                    order = [1,0,2]
                    wd_ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], \
                                   loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                                   ncol = 3,frameon = False,columnspacing = plt_param.LEGcol_space);    
                    ###############################################################
                    wd_ax.set_ylim(WDmin_ext,WDmax_ext); wd_ax.set_yticks(WDtick_ext);
                    wd_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                    ################################################################
                    wd_ax.set_ylabel('Wind Direction, deg ${}^{\circ}$',{'fontsize': plt_param.LABELfs});
                    wd_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                    ###############################################################
                    #insert data date as title
                    T = wd_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                    T.set_position(plt_param.Tloc);
                    
                    """ ###########################################################
                    # Plot the Model Wind Direction Biases
                    ########################################################### """
                    # initate 0/0 line
                    bias_ax.plot([Xmin_ext,Xmax_ext],[0,0],color='black',linewidth=1.5)
                    ###############################################################
                    # Define the Bias (Observations - Measurements)
                    WDbias,RMSEwd,BIASwd = wd_bias(MEASwd_ext_mu,MDLwd_ext)
                    bias_ax.plot(MDLtime_ext,WDbias,'^-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Mdl - Obs')
                    bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                    ###############################################################
                    bias_ax.legend(loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                                   ncol = 1,frameon = False);      
                    ###############################################################
                    #bias_ax.set_xticks(Xtick_ext); bias_ax.set_xticklabels(Xtick_ext_lab,{'fontsize': plt_param.TICKfs});
                    bias_ax.set_xticks(Xtick_ext); bias_ax.set_xticklabels(Xtick_ext_lab,fontsize=plt_param.TICKfs);
                    bias_ax.set_xlim(Xmin_ext,Xmax_ext)
                    ###############################################################
                    bias_ax.set_xlabel('Hour, UTC',{'fontsize': plt_param.LABELfs});
                    bias_ax.xaxis.set_label_coords(plt_param.Xlab_pad_dual[0],plt_param.Xlab_pad_dual[1]);
                    ###############################################################
                    #Define Bounds of Wind Speed - Model Bias to Determine Y-Bounds
                    WDb_min,WDb_max,WDb_tick = wd_blim((WDbias));
                    bias_ax.set_ylim(WDb_min,WDb_max); bias_ax.set_yticks(WDb_tick);
                    bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                    ###############################################################
                    bias_ax.set_ylabel('Model Error',{'fontsize': plt_param.LABELfs});    
                    bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                    ###############################################################
                    #insert data statistics as title
                    BT = bias_ax.set_title('RMSE: ' + str(RMSEwd)+' deg '+r'$\mathrm{}^{\circ}$ | Mean Bias: '+\
                                      str(BIASwd)+' deg '+r'$\mathrm{}^{\circ}$',\
                                      {'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                    BT.set_position(plt_param.Bloc_left);
                    ###############################################################
                    #-> annotate measurement heights
                    wd_ax.annotate(str(MEAShubhgt_value)+' m AGL',xy=(0.13,0.895), xycoords='figure fraction',color='black',\
                                     fontsize=plt_param.LEGfs,weight='bold')
                    ###############################################################
                    wd_fig.subplots_adjust(left = 0.125,right = 0.935,bottom = 0.155, top = 0.925)
                    ###############################################################
                    bias_ax.annotate(Date_start_label_ext,xy = plt_param.metDATEpad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                                #fontsize = plt_param.IDfs,bbox = dict(boxstyle='square',fc='w'))
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
                    imax = wd_fig.add_axes(plt_param.metLOGOaxcoords)
                    imax.set_axis_off()
                    imax.imshow(logo, aspect="equal")
                    ###################################################################
                    #Save the Developed Figure
                    FILout = os.path.join(IMGpth_ext,FILpfx+MDLname+'_WDhub_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                    wd_fig.savefig(FILout,dpi=300); plt.close('all')
    
                    """ Extended Power Forecast """
                    pow_fig, (pow_ax,bias_ax) = plt.subplots(2,1,sharex = True, figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                         facecolor = 'w',edgecolor ='w',clear='True')
                    ###############################################################
                    # Plot the Power Time Histories
                    #####################################
                    #Observations 
                    pow_ax.plot(MDLtime_ext,MEASpow_ext_mu,'s-',color='gray',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Obs ('+AVGlabel+' Avg)')
                    pow_ax.plot(MEAStime,MEASpow,'.',color='black',markersize=1,linewidth=0,label='Obs')
                    pow_ax.plot(MDLtime_ext,MDLpow_ext,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname+' ('+str(MDLhgt_ref_value)+' m)')
                    pow_ax.add_collection(LineCollection(MEASpow_ext_up_segs,colors=(plt_param.POWramp_up_color),linewidths=(plt_param.PLTlwidth),label='Up ramps'))
                    pow_ax.add_collection(LineCollection(MEASpow_ext_down_segs,colors=(plt_param.POWramp_down_color),linewidths=(plt_param.PLTlwidth),label='Down ramps'))
                    pow_ax.add_collection(LineCollection(MDLpow_ext_up_segs,colors=(plt_param.POWramp_up_color),linewidths=(plt_param.PLTlwidth),label=' '))
                    pow_ax.add_collection(LineCollection(MDLpow_ext_down_segs,colors=(plt_param.POWramp_down_color),linewidths=(plt_param.PLTlwidth),label=' '))
                    pow_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                    ###################################################################
                    handles, labels = pow_ax.get_legend_handles_labels()
                    # override line label 5 to be a dummy white line
                    handles[5] = Line2D([0],[0],color="w")
                    order = [1,5,0,3,2,4]
                    pow_ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], \
                                   loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                                   ncol = 3,frameon = False,columnspacing = plt_param.LEGcol_space);    
                    ###############################################################
                    pow_ax.set_ylim(POWmin_ext,POWmax_ext); pow_ax.set_yticks(POWtick_ext);
                    pow_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                    ###############################################################
                    pow_ax.set_ylabel('Wind Capacity Factor',{'fontsize': plt_param.LABELfs});
                    pow_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                    ###############################################################
                    #insert station info as title
                    T = pow_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                    T.set_position(plt_param.Tloc);
                                 
                    """ ###########################################################
                    # Plot the Model Power Biases
                    ########################################################### """
                    # initate 0/0 line
                    bias_ax.plot([Xmin_ext,Xmax_ext],[0,0],color='black',linewidth=1.5)
                    ###############################################################
                    # Define the Bias (Observations - Measurements)
                    POWbias = MDLpow_ext - MEASpow_ext_mu
                    bias_ax.plot(MDLtime_ext,POWbias,'^-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Mdl - Obs')
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
                    #Define Bounds of Power - Model Bias to Determine Y-Bounds
                    #POWb_min,POWb_max,POWb_tick = atmo_bound_lim((POWbias),plt_param.POWbias_bse,plt_param.POWbmin,\
                    #                                 plt_param.POWbmax,plt_param.POWbias_bse)
                    POWb_min=-1.1;POWb_max=1.1;POWb_tick=np.arange(-1,1.1,.5);
                    bias_ax.set_ylim(POWb_min,POWb_max); bias_ax.set_yticks(POWb_tick);
                    bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                    ###############################################################
                    bias_ax.set_ylabel('Model Error',{'fontsize': plt_param.LABELfs});    
                    bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                    ###############################################################
                    #input basic statistics
                    RMSEpow,BIASpow = atmo_bias(POWbias)
                    #insert data statistics as title
                    BT = bias_ax.set_title('RMSE: ' + str(RMSEpow)+' | Mean Bias: '+\
                                          str(BIASpow)+' ',\
                                          {'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5) 
                    BT.set_position(plt_param.Bloc_left);
                    ###############################################################
                    #-> annotate measurement heights
                    pow_ax.annotate(str(MEAShubhgt_value)+' m AGL',xy=(0.13,0.895), xycoords='figure fraction',color='black',\
                                     fontsize=plt_param.LEGfs,weight='bold')
                    ###############################################################
                    pow_fig.subplots_adjust(left = 0.125,right = 0.935,bottom = 0.155, top = 0.925)
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
                    imax = pow_fig.add_axes(plt_param.metLOGOaxcoords)
                    imax.set_axis_off()
                    imax.imshow(logo, aspect="equal")
                    #######################################################################
                    #Save the Developed Figure
                    FILout = os.path.join(IMGpth_ext,FILpfx+MDLname+'_POW_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                    pow_fig.savefig(FILout,dpi=300); plt.close('all');
                    ###############################################################
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
                Xtick = np.arange(MDLtime[0],MDLtime[0]+(FXlen - 1)*HRsec+1e-5,plt_param.Xdel)
                Xtick_lab = [time.strftime('%H', time.gmtime(time_ref)) for time_ref in Xtick]
                Time_range_dates = [basedate+timedelta(seconds=(time_ref)) for time_ref in Xtick]
                Date_start_label=Time_range_dates[0].strftime('%Y-%m-%d')
                Date_end_label=Time_range_dates[-1].strftime('%Y-%m-%d')
                    
                """ Perform Model Averaging for Plotting and Statistical Analyses """
                MEASws_mu = atmo_meas2hr_mid(MEASws,None,MEAStime,MDLtime,AVGmins=AVGper,func='median');
                MEASpow_mu = atmo_meas2hr_mid(MEASpow,None,MEAStime,MDLtime,AVGmins=AVGper,func='median');
                MEASwd_mu = wd_meas2hr_mid(MEASwd,None,MEAStime,MDLtime,AVGmins=AVGper,func='median');
                MEASws_up_segs, MEASws_down_segs = ramp_up_down_segs(MDLtime,MEASws_mu,threshold=plt_param.WSthreshold)
                MEASpow_up_segs, MEASpow_down_segs = ramp_up_down_segs(MDLtime,MEASpow_mu,threshold=plt_param.POWthreshold)
                
                """ Define Relevant Plotting Specifics """
                # Define Proper Indices to Determine Axes Parameters
                meas_oi = np.logical_and(MEAStime >= MDLtime[0], MEAStime <= MDLtime[-1]) 
                ###################################################################
                WSmin,WSmax,WStick = atmo_bound_lim((MEASws[meas_oi]),\
                                               plt_param.WSrnd_bse,plt_param.WSmin,plt_param.WSmax,\
                                               plt_param.WSrnd_bse)
                #POWmin,POWmax,POWtick = atmo_bound_lim((MEASpow[meas_oi]),\
                #                               plt_param.POWrnd_bse,plt_param.POWmin,plt_param.POWmax,\
                #                               plt_param.POWrnd_bse)
                POWmin=plt_param.POWmin;POWmax=plt_param.POWmax;POWtick=np.arange(0,plt_param.POWmax,plt_param.POWrnd_bse);
                WDmin,WDmax,WDtick = wd_ylim((MEASwd[meas_oi]));
    
                """ ###############################################################
                Wind Speed Plots (instruct to share x-axis)
                ############################################################### """
                ws_fig, (ws_ax,bias_ax) = plt.subplots(2,1,sharex = True, figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                         facecolor = 'w',edgecolor ='w',clear='True')
                ###################################################################
                # Plot the Wind Speed Time Histories
                ####################################
                #Observations -- Only Plotting the QC'd Values
                ws_ax.plot(MDLtime,MEASws_mu,'s-',color='gray',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Obs ('+AVGlabel+' Avg)')
                ws_ax.plot(MEAStime,MEASws,'.',color='black',markersize=1,linewidth=0,label='Obs')
                ws_ax.plot(MDLtime,MDLtime*np.nan,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname)
                ws_ax.add_collection(LineCollection(MEASws_up_segs,colors=(plt_param.WSramp_up_color),linewidths=(plt_param.PLTlwidth),label='Up ramps'))
                ws_ax.add_collection(LineCollection(MEASws_down_segs,colors=(plt_param.WSramp_down_color),linewidths=(plt_param.PLTlwidth),label='Down ramps'))
                ws_ax.plot(MDLtime,MDLtime*np.nan,'-',color=plt_param.WSramp_up_color,markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=' ')
                ws_ax.plot(MDLtime,MDLtime*np.nan,'-',color=plt_param.WSramp_down_color,markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='')
                ws_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                handles, labels = ws_ax.get_legend_handles_labels()
                # override line label 5 to be a dummy white line
                handles[5] = Line2D([0],[0],color="w")
                order = [1,5,0,3,2,4]
                ws_ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], \
                               loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                               ncol = 3,frameon = False,columnspacing = plt_param.LEGcol_space);    
                ###################################################################
                ws_ax.set_ylim(WSmin,WSmax); ws_ax.set_yticks(WStick);
                ws_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                ws_ax.set_ylabel('Wind Speed, m $\mathrm{s}^{-1}$',{'fontsize': plt_param.LABELfs});
                ws_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                ###################################################################
                #insert data date as title
                T = ws_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                T.set_position(plt_param.Tloc);
                        
                """ ###############################################################
                Define and Plot the Model Wind Speed Biases
                ############################################################### """
                # initate 0/0 line
                bias_ax.plot([Xmin,Xmax],[0,0],color='black',linewidth=1.5)
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
                bias_ax.text(np.nanmean([Xmin,Xmax]),0,MDLname+' Data Not Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
                ###################################################################
                #Define Bounds of Wind Speed - Model Bias to Determine Y-Bounds
                bias_ax.set_ylim(-3,3); bias_ax.set_yticks(np.arange(-3,3.1,1));
                bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                bias_ax.set_ylabel('Model Error',{'fontsize': plt_param.LABELfs});    
                bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                ###################################################################
                #insert data statistics as title
                BT = bias_ax.set_title('RMSE: N/A',{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                BT.set_position(plt_param.Bloc_left);
                ###################################################################
                #-> annotate measurement heights
                ws_ax.annotate(str(MEAShubhgt_value)+' m AGL',xy=(0.13,0.895), xycoords='figure fraction',color='black',\
                                 fontsize=plt_param.LEGfs,weight='bold')
                ###################################################################
                ws_fig.subplots_adjust(left = 0.125,right = 0.935,bottom = 0.155, top = 0.925)
                ###################################################################
                bias_ax.annotate(Date_start_label,xy = plt_param.metDATEpad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                                #fontsize = plt_param.IDfs,bbox = dict(boxstyle='square',fc='w'))
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
                imax = ws_fig.add_axes(plt_param.metLOGOaxcoords)
                imax.set_axis_off()
                imax.imshow(logo, aspect="equal")
                ###################################################################
                #Save the Developed Figure
                FILout = os.path.join(IMGpth_model,FILpfx+MDLname+'_WShub_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                ws_fig.savefig(FILout,dpi=300); plt.close('all');
                        
                """ ###############################################################
                Wind Direction Plots (instruct to share x-axis)
                ############################################################### """
                wd_fig, (wd_ax,bias_ax) = plt.subplots(2,1,sharex = True, figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                         facecolor = 'w',edgecolor ='w',clear='True')
                ###################################################################
                # Plot the Wind Direction Time Histories
                ###################################################################
                wd_ax = wd_plt(MDLtime,MEASwd_mu,wd_ax,'s-','gray','Obs ('+AVGlabel+' Avg)',marker_size=plt_param.PLTmsize)
                wd_ax = wd_plt(MEAStime,MEASwd,wd_ax,'.','black','Obs',marker_size=1,line_width=0)
                wd_ax.plot(MDLtime,MDLtime*np.nan,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname)
                ###################################################################
                wd_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                handles, labels = wd_ax.get_legend_handles_labels()
                order = [1,0,2]
                wd_ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], \
                                   loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                                   ncol = 3,frameon = False,columnspacing = plt_param.LEGcol_space);    
                ###################################################################
                wd_ax.set_ylim(WDmin,WDmax); wd_ax.set_yticks(WDtick);
                wd_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                wd_ax.set_ylabel('Wind Direction, deg ${}^{\circ}$',{'fontsize': plt_param.LABELfs});
                wd_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                ###################################################################
                #insert data date as title
                T = wd_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                T.set_position(plt_param.Tloc);
                        
                """ ###############################################################
                # Plot the Model Wind Direction Biases
                ############################################################### """
                # initate 0/0 line
                bias_ax.plot([Xmin,Xmax],[0,0],color='black',linewidth=1.5)
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
                bias_ax.text(np.nanmean([Xmin,Xmax]),0,MDLname+' Data Not Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
                ###################################################################
                #Define Bounds of Wind Speed - Model Bias to Determine Y-Bounds
                bias_ax.set_ylim(-3,3); bias_ax.set_yticks(np.arange(-3,3.1,1));
                bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                bias_ax.set_ylabel('Model Error',{'fontsize': plt_param.LABELfs});    
                bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                ###################################################################
                #insert data statistics as title
                BT = bias_ax.set_title('RMSE: N/A',{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                BT.set_position(plt_param.Bloc_left);
                ###################################################################
                #-> annotate measurement heights
                wd_ax.annotate(str(MEAShubhgt_value)+' m AGL',xy=(0.13,0.895), xycoords='figure fraction',color='black',\
                                 fontsize=plt_param.LEGfs,weight='bold')
                ###################################################################
                wd_fig.subplots_adjust(left = 0.125,right = 0.935,bottom = 0.155, top = 0.925)
                ###################################################################
                bias_ax.annotate(Date_start_label,xy = plt_param.metDATEpad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                                #fontsize = plt_param.IDfs,bbox = dict(boxstyle='square',fc='w'))
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
                imax = wd_fig.add_axes(plt_param.metLOGOaxcoords)
                imax.set_axis_off()
                imax.imshow(logo, aspect="equal")
                ###################################################################
                #Save the Developed Figure
                FILout = os.path.join(IMGpth_model,FILpfx+MDLname+'_WDhub_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                wd_fig.savefig(FILout,dpi=300); plt.close('all')
                ###################################################################
                """ ###############################################################
                Power Profile (instruct to share x-axis)
                ############################################################### """
                pow_fig, (pow_ax,bias_ax) = plt.subplots(2,1,sharex = True, figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                         facecolor = 'w',edgecolor ='w',clear='True')
                ####################################################################
                # Plot the Power Time Histories
                #####################################
                #Observations 
                pow_ax.plot(MDLtime,MEASpow_mu,'s-',color='gray',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Obs ('+AVGlabel+' Avg)')
                pow_ax.plot(MEAStime,MEASpow,'.',color='black',markersize=1,linewidth=0,label='Obs')
                pow_ax.plot(MDLtime,MDLtime*np.nan,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname)
                pow_ax.add_collection(LineCollection(MEASpow_up_segs,colors=(plt_param.POWramp_up_color),linewidths=(plt_param.PLTlwidth),label='Up ramps'))
                pow_ax.add_collection(LineCollection(MEASpow_down_segs,colors=(plt_param.POWramp_down_color),linewidths=(plt_param.PLTlwidth),label='Down ramps'))
                pow_ax.plot(MDLtime,MDLtime*np.nan,'-',color=plt_param.WSramp_up_color,markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=' ')
                pow_ax.plot(MDLtime,MDLtime*np.nan,'-',color=plt_param.WSramp_down_color,markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='')
                pow_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                handles, labels = pow_ax.get_legend_handles_labels()
                # override line label 5 to be a dummy white line
                handles[5] = Line2D([0],[0],color="w")
                order = [1,5,0,3,2,4]
                pow_ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], \
                               loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                               ncol = 3,frameon = False,columnspacing = plt_param.LEGcol_space);    
                ###################################################################
                pow_ax.set_ylim(POWmin,POWmax); pow_ax.set_yticks(POWtick);
                pow_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                pow_ax.set_ylabel('Wind Capacity Factor',{'fontsize': plt_param.LABELfs});
                pow_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                ###################################################################
                #insert data date as title
                T = pow_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                T.set_position(plt_param.Tloc);
                             
                """ ###############################################################
                # Plot the Model Power Biases
                ############################################################### """
                # initate 0/0 line
                bias_ax.plot([Xmin,Xmax],[0,0],color='black',linewidth=1.5)
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
                bias_ax.text(np.nanmean([Xmin,Xmax]),0,MDLname+' Data Not Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
                ###################################################################
                #Define Bounds of Wind Speed - Model Bias to Determine Y-Bounds
                bias_ax.set_ylim(-1.1,1.1); bias_ax.set_yticks(np.arange(-1,1.1,.5));
                bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                bias_ax.set_ylabel('Model Error',{'fontsize': plt_param.LABELfs});    
                bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                ###################################################################
                #insert data statistics as title
                BT = bias_ax.set_title('RMSE: N/A',{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                BT.set_position(plt_param.Bloc_left);
                ###################################################################
                #-> annotate measurement heights
                pow_ax.annotate(str(MEAShubhgt_value)+' m AGL',xy=(0.13,0.895), xycoords='figure fraction',color='black',\
                                 fontsize=plt_param.LEGfs,weight='bold')
                ###################################################################
                pow_fig.subplots_adjust(left = 0.125,right = 0.935,bottom = 0.155, top = 0.925)
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
                imax = pow_fig.add_axes(plt_param.metLOGOaxcoords)
                imax.set_axis_off()
                imax.imshow(logo, aspect="equal")
                #######################################################################
                #Save the Developed Figure
                FILout = os.path.join(IMGpth_model,FILpfx+MDLname+'_POW_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                pow_fig.savefig(FILout,dpi=300); plt.close('all');
