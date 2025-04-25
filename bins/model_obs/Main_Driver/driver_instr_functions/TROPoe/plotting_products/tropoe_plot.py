#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The Purpose of this Script is to Output the TROPoe Products from the PSL website

###############################################################################
Created on Nov 25 13:31:49 2023
@author: dmurray (after jduncan)
###############################################################################
Outputted Plots Follow:
instrument_whatisplotted_Sid.date.png
**Imperative that the "." is the first period -- this is used to strip the **
**date from the file so that duplicate images are not run **
###############################################################################
Aside: Differences Between Pcolor() and Pcolormesh()
    While pcolor returns a PolyCollection, pcolormesh returns a QuadMesh. 
    The latter is more specialized for the given purpose and thus is faster. 
    It should almost always be preferred. -- plt.pcolor(X,Y,TEMP)
    ###########################################################################
    - plotting methods include:
    t_img = t_ax.pcolormesh(X,Y,TEMP,vmin=Vmin,vmax=Vmax,shading='gouraud')   
    t_img = t_ax.imshow(TEMP,extent=IMext,vmin=Vmin,vmax=Vmax,interpolation='bilinear',\
                          origin='lower',aspect='auto')                    

Using MatPlotLib's Pcolormesh to Visualize the TROPoe Fields.
- Would like to use some type of interpolation field (i.e. Gouraud) to improve
  optics. However, when applying Gouraud (using pcolormesh) the areas 
 neighbouring NaNs appear with some weird grey structure indicating that the 
  grid box could not be properly implemented
- This describes the issue:
    https://github.com/matplotlib/matplotlib/issues/8802

###############################################################################
# No QC Information was Provided for these Datastreams
        
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as image
import matplotlib.colors as colors
import os,time,warnings
from datetime import datetime,timedelta
#####################################
from coord_distance import dist2coord
#####################################
from atmo_meas_to_grid import atmo_meas2grid
from atmo_lin_interp2D import atmo_lininterp2_lev2D
######################################
from time_tuple_cat import time_cat
from atmo_temperature import calc_theta
from atmo_temperature import degC_to_K
from atmo_spatial_tuple_cat import atmo_spatial_cat
from atmo_tuple_cat import atmo_cat
from atmo_spatial_correction import atmo_sptl_correct
###################################################
from atmo_ylim_plot import atmo_ylim
from atmo_spatial_ylim import atmo_sptl_ylim
from atmo_bias_spatial_ylim import atmo_sptlbias_ylim
from avgmins_to_label import avgmins_to_label
#####################################################
from cmaps import get_var_cmap
from label_funcs import gridpoint_to_label
from format_funcs import format_lat_lon_alt
from matplotlib.colors import from_levels_and_colors
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
AVGper = getattr(plt_param,'TROPoe_avgmins',60)
ccolor = plt_param.PLTccolor

def tropoe_plot(MEASin,MDLin,OUTinfo,version='ASSIST'):
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
    tropoe_plot(TROPoe_plt,HRRRw_plt,FILout)
    ###########################################################################
    MDLin ::  HRRRv4_plt = (HRRRv4_ini,HRRRv4_loc,HRRRv4_xy,HRRRv4_trp)
    OUTinfo :: FILout = (IMGout,MDLname,Sid,doi)
    ####################################################################### """
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
    """ Define Measurement Location and Model Grid Point Structure """
    ###########################################################################
    # measurement data not defined -- load in location from station parameter file
    MEASlon = getattr(station_param,Sid+'lon')
    MEASlat = getattr(station_param,Sid+'lat')
    MEASalt = getattr(station_param,Sid+'alt')
    MEASid = getattr(station_param,Sid+'name')
    MEASprojid = getattr(station_param,Sid+'projID',MEASid)
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

    # Declare Measurement Heights
    MEAShgt = MEASin[1][1][1];  # converted to meters in tropoe_main_ingest
    #MAXhgt = MEAShgt[-1]
    MAXhgt = 3000
    # Define Default cloud base height (not from ceilometer)
    if DATE < "20240112":
        DEFcbh=2000
    else:
        DEFcbh=1000
    
    # Define Relevant Indices for Reference in Concatenation #
    T_cat = 0; RH_cat = 1; 
    TD_cat = 2; MX_cat = 3;
    THETA_cat = 0;
    CBHind_cat = 0; GAMMAind_cat = 1; RMSRind_cat=2; LWPind_cat=3;
    ###########################################################################
    # For Spatial Concatenation -- Define the Time Axis within the Input Tuple
    Tax = 1;
        
    #-> Combine Relevant Measurement Tuples 
    MEAStime = time_cat(MEASin[1][0],MEASin[1][1],MEASin[1][2],MEASin[1][3])
    ###########################################################################
    MEAStrh_in = MEASin[2]; 
    MEAStheta_in = MEASin[3]; 
    MEASmisc_in = MEASin[4]; 
    #####################
    MEASt = atmo_spatial_cat(MEAStrh_in[0],MEAStrh_in[1],MEAStrh_in[2],MEAStrh_in[3],T_cat,Tax)
    MEAStd = atmo_spatial_cat(MEAStrh_in[0],MEAStrh_in[1],MEAStrh_in[2],MEAStrh_in[3],TD_cat,Tax)
    MEASrh = atmo_spatial_cat(MEAStrh_in[0],MEAStrh_in[1],MEAStrh_in[2],MEAStrh_in[3],RH_cat,Tax)
    MEASmx = atmo_spatial_cat(MEAStrh_in[0],MEAStrh_in[1],MEAStrh_in[2],MEAStrh_in[3],MX_cat,Tax)
    MEAStheta = atmo_spatial_cat(MEAStheta_in[0],MEAStheta_in[1],MEAStheta_in[2],MEAStheta_in[3],THETA_cat,Tax)
    MEAScbh = atmo_cat(MEASmisc_in[0],MEASmisc_in[1],MEASmisc_in[2],MEASmisc_in[3],CBHind_cat)*1000 # to meters
    MEASgamma = atmo_cat(MEASmisc_in[0],MEASmisc_in[1],MEASmisc_in[2],MEASmisc_in[3],GAMMAind_cat)
    MEASrmsr = atmo_cat(MEASmisc_in[0],MEASmisc_in[1],MEASmisc_in[2],MEASmisc_in[3],RMSRind_cat)
    MEASlwp = atmo_cat(MEASmisc_in[0],MEASmisc_in[1],MEASmisc_in[2],MEASmisc_in[3],LWPind_cat)
    
    # set values to nan when certain conditions are met
    for t in range(len(MEASgamma)):
       if not np.logical_and(MEASgamma[t] == 1, MEASrmsr[t] < 5):
           MEASt[:,t] = np.nan
           MEASmx[:,t] = np.nan
           MEAStheta[:,t] = np.nan

    # set the cloud base height to nan when lwp is < 5
    MEAScbh[np.where(MEASlwp < 5)] = np.nan;

    #  For ASSIST plots, create a version of the plotted data with values above cbh set to NaN
    if (version == 'ASSIST'):
        MEASt_nan = np.copy(MEASt); MEASmx_nan = np.copy(MEASmx); MEAStheta_nan = np.copy(MEAStheta);
        for t in range(len(MEAStime)):
            #only when LWP > 5
            if MEASlwp[t]>=5:
                #idx=np.where((MEAShgt>MEAScbh[t]) | ((MEAShgt>0) & (MEAScbh[t]==DEFcbh)))[0]
                idx=np.where(MEAShgt>MEAScbh[t])[0]
                MEASt_nan[idx,t]=np.nan
                MEASmx_nan[idx,t]=np.nan
                MEAStheta_nan[idx,t]=np.nan
                idx2 =np.where((MEAShgt>0) & (MEAScbh[t]==DEFcbh))[0]
                MEASt_nan[idx2,t]=np.nan
                MEASmx_nan[idx2,t]=np.nan
                MEAStheta_nan[idx2,t]=np.nan

    ###################################################################
    # Perform Correction to Create Equally Spaced Data Arrays #
    ###################################################################
    MEAStme_base, MEAShgt_base = np.meshgrid(MEAStime,MEAShgt) #define reference meshgrid
    MEAStme_mesh, MEAShgt_mesh, MEASt = atmo_sptl_correct(MEAStme_base,MEAShgt_base,MEASt);
    _, _, MEAStd = atmo_sptl_correct(MEAStme_base,MEAShgt_base,MEAStd);
    _, _, MEASrh = atmo_sptl_correct(MEAStme_base,MEAShgt_base,MEASrh);
    _, _, MEASmx = atmo_sptl_correct(MEAStme_base,MEAShgt_base,MEASmx);
    _, _, MEAStheta = atmo_sptl_correct(MEAStme_base,MEAShgt_base,MEAStheta);
    _, _, MEAScbh = atmo_sptl_correct(MEAStme_base,MEAShgt_base,MEAScbh);
    if (version == 'ASSIST'):
        _, _, MEASt_nan = atmo_sptl_correct(MEAStme_base,MEAShgt_base,MEASt_nan);
        _, _, MEASmx_nan = atmo_sptl_correct(MEAStme_base,MEAShgt_base,MEASmx_nan);
        _, _, MEAStheta_nan = atmo_sptl_correct(MEAStme_base,MEAShgt_base,MEAStheta_nan);

    # create versions of cbh for plotting based on whether it is the default or from ceilometer
    MEAScbh_def = np.copy(MEAScbh)
    MEAScbh_def[np.where(MEAScbh != DEFcbh)] = np.nan;
    MEAScbh_nondef = np.copy(MEAScbh)
    MEAScbh_nondef[np.where(MEAScbh==DEFcbh)] = np.nan;
    
    """#######################################################################
    Define some stuff for the plots
    ########################################################################## """
    # set up the colormaps
    cmap_temp = get_var_cmap('temp')
    cmap_theta = get_var_cmap('temp')
    cmap_mx = get_var_cmap('mx')
    cmap_bias_temp = get_var_cmap('bias_temp')
    cmap_bias_theta = get_var_cmap('bias_temp')
    cmap_bias_mx = get_var_cmap('bias_moisture')

    # set up the variable scale ranges
    TEMPmin=plt_param.TROPoe_TEMPmin;TEMPmax=plt_param.TROPoe_TEMPmax;TEMPtick=plt_param.TROPoe_TEMPtick;
    TEMPclabels=np.arange(TEMPmin,TEMPmax,plt_param.TEMPcint)
    TEMPmin_ext=TEMPmin;TEMPmax_ext=TEMPmax;TEMPtick_ext=TEMPtick;TEMPclabels_ext=TEMPclabels;
    TEMPb_min=plt_param.TROPoe_TEMPbmin;TEMPb_max=plt_param.TROPoe_TEMPbmax;TEMPb_tick=plt_param.TROPoe_TEMPbtick;
    TEMPb_bounds=np.arange(plt_param.TROPoe_TEMPbmin,plt_param.TROPoe_TEMPbmax+1,2);
    TEMPb_tick=np.arange(-8,9,4)

    MIXRmin=plt_param.TROPoe_MIXRmin;MIXRmax=plt_param.TROPoe_MIXRmax;MIXRtick=plt_param.TROPoe_MIXRtick;
    MIXRclabels=np.arange(MIXRmin,MIXRmax,plt_param.MIXRcint)
    MIXRmin_ext=MIXRmin;MIXRmax_ext=MIXRmax;MIXRtick_ext=MIXRtick;MIXRclabels_ext=MIXRclabels;
    MIXRb_min=plt_param.TROPoe_MIXRbmin;MIXRb_max=plt_param.TROPoe_MIXRbmax;MIXRb_tick=plt_param.TROPoe_MIXRbtick;
    MIXRb_bounds=np.arange(plt_param.TROPoe_MIXRbmin,plt_param.TROPoe_MIXRbmax+1,1);
    MIXRb_tick=np.arange(-6,6,2);

    THETAmin=plt_param.TROPoe_THETAmin;THETAmax=plt_param.TROPoe_THETAmax;THETAtick=plt_param.TROPoe_THETAtick;
    THETAclabels=np.arange(THETAmin,THETAmax,plt_param.THETAcint)
    THETAmin_ext=THETAmin;THETAmax_ext=THETAmax;THETAtick_ext=THETAtick;THETAclabels_ext=THETAclabels;
    THETAb_min=plt_param.TROPoe_THETAbmin;THETAb_max=plt_param.TROPoe_THETAbmax;THETAb_tick=plt_param.TROPoe_THETAbtick;
    THETAb_bounds=np.arange(plt_param.TROPoe_THETAbmin,plt_param.TROPoe_THETAbmax+1,1);
    THETAb_tick=np.arange(-8,9,4);

    GMLcbhlabel = 'Cloud Base Height (from GML ceilometer)'
    DEFcbhlabel = 'Cloud Base Height (default value, true value unknown)'

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
                ####################################################
                MDLt_ini = MDLtpmx[INIref][INTPoi]
                ###################################################################
                # gridded interpolated tpmx information (Tint,Tsfc_int,Pint,Psfc_int,MXrat_int,MXrat_sfc_int,Prate_int)
                MDLt = MDLt_ini[0]; MDLp = MDLt_ini[2]; MDLmx = MDLt_ini[4];
                # data is currently shaped as: (forecast_hour,sites,hgt)
                ###################################################################
                # extract from relevant grid point
                MDLt = MDLt[:,Gref,:]; MDLp = MDLp[:,Gref,:]; MDLmx = MDLmx[:,Gref,:];
                MDLtheta = calc_theta(degC_to_K(MDLt),MDLp)
                ###################################################################
                # Account for the way python interprets matrices (i.e. row x column)
                # ---> transpose inputted model data
                MDLt = np.transpose(MDLt); MDLmx = np.transpose(MDLmx);
                MDLtheta = np.transpose(MDLtheta);
                # now correctly (hgt x time) and consistent with the measured data
                ###################################################################
                # Interpolate the MDL vars to the MEAS heights for comparison
                MDLMEASt = atmo_lininterp2_lev2D(MDLt,MDLhgt,MEAShgt)
                MDLMEASmx = atmo_lininterp2_lev2D(MDLmx,MDLhgt,MEAShgt)
                MDLMEAStheta = atmo_lininterp2_lev2D(MDLtheta,MDLhgt,MEAShgt)
                ###################################################################
                """ Define Relevant Plotting Specifics 
                --> Specifics will vary depending on whether an extended forecast 
                --> is being produced
                ###################################################################
                """
                ###################################################################
                #MEAStme_base, MEAShgt_base = np.meshgrid(MEAStime,MEAShgt) #define reference meshgrid
                ###################################################################
                # Perform Correction to Create Equally Spaced Data Arrays #
                ###################################################################
                #MEAStme_mesh = MEAStme_base; MEAShgt_mesh = MEAShgt_base
                #MEAStme_mesh, MEAShgt_mesh, MEASt = atmo_sptl_correct(MEAStme_base,MEAShgt_base,MEASt);
                #_, _, MEAStd = atmo_sptl_correct(MEAStme_base,MEAShgt_base,MEAStd);
                #_, _, MEASrh = atmo_sptl_correct(MEAStme_base,MEAShgt_base,MEASrh);
                #_, _, MEASmx = atmo_sptl_correct(MEAStme_base,MEAShgt_base,MEASmx);
                #_, _, MEAStheta = atmo_sptl_correct(MEAStme_base,MEAShgt_base,MEAStheta);
                #_, _, MEAScbh = atmo_sptl_correct(MEAStme_base,MEAShgt_base,MEAScbh);
                #MEAScbh_def = np.copy(MEAScbh)
                #MEAScbh_def[np.where(MEAScbh != DEFcbh)] = np.nan;
                #MEAScbh_nondef = np.copy(MEAScbh)
                #MEAScbh_nondef[np.where(MEAScbh==DEFcbh)] = np.nan;
                #if (version == 'ASSIST'):
                #    _, _, MEASt_nan = atmo_sptl_correct(MEAStme_base,MEAShgt_base,MEASt_nan);
                #    _, _, MEASmx_nan = atmo_sptl_correct(MEAStme_base,MEAShgt_base,MEASmx_nan);
                #    _, _, MEAStheta_nan = atmo_sptl_correct(MEAStme_base,MEAShgt_base,MEAStheta_nan);
    
                ###################################################################
                # maximum and minimum measurement height 
                HGTmin,HGTmax,HGTtick = atmo_ylim((np.asarray([0,MAXhgt])),plt_param.SHGTrnd_bse,plt_param.TICKint);
                HGTtick_label = [str(ht/1000) for ht in HGTtick]
                ###################################################################
                # Define the Max Height Ind
                MEAShgt_max = np.where(MEAShgt > MAXhgt)[0]
                MDLhgt_max = np.where(MDLhgt > MAXhgt)[0]
    
                ###################################################################
                if len(MDLtime) > FXlen:
                    """ Define Bounds for Extended Forecast Period """
                    ###############################################################
                    # Adjust Model Variables #
                    MDLt_ext = MDLt; MDLmx_ext = MDLmx; MDLtheta_ext = MDLtheta;
                    MDLMEASt_ext = MDLMEASt; MDLMEASmx_ext = MDLMEASmx;
                    MDLMEAStheta_ext = MDLMEAStheta;
                    MDLtime_ext = MDLtime;
                    ###############################################################
                    """ Perform Model Averaging for Plotting and Statistical Analyses """
                    MDLtme_ext_mesh, MDLhgt_ext_mesh = np.meshgrid(MDLtime_ext,MDLhgt); #define reference meshgrid
                    MDLMEAStme_ext_mesh, MDLMEAShgt_ext_mesh = np.meshgrid(MDLtime_ext,MEAShgt)
                    ###############################################################
                    if (version == 'ASSIST'):
                        MEASt_ext_mu = atmo_meas2grid(MEASt_nan,(MEAShgt_mesh,MEAStme_mesh),(MEAShgt,MDLtime_ext),AVGmins=AVGper)
                        MEASmx_ext_mu = atmo_meas2grid(MEASmx_nan,(MEAShgt_mesh,MEAStme_mesh),(MEAShgt,MDLtime_ext),AVGmins=AVGper)
                        MEAStheta_ext_mu = atmo_meas2grid(MEAStheta_nan,(MEAShgt_mesh,MEAStme_mesh),(MEAShgt,MDLtime_ext),AVGmins=AVGper)
                    else:
                        MEASt_ext_mu = atmo_meas2grid(MEASt,(MEAShgt_mesh,MEAStme_mesh),(MEAShgt,MDLtime_ext),AVGmins=AVGper)
                        MEASmx_ext_mu = atmo_meas2grid(MEASmx,(MEAShgt_mesh,MEAStme_mesh),(MEAShgt,MDLtime_ext),AVGmins=AVGper)
                        MEAStheta_ext_mu = atmo_meas2grid(MEAStheta,(MEAShgt_mesh,MEAStme_mesh),(MEAShgt,MDLtime_ext),AVGmins=AVGper)
                    ###############################################################
                    ext_meas_oi = np.logical_and(MEAShgt_mesh <= MAXhgt,np.logical_and(MEAStme_mesh >= MDLtime_ext[0], MEAStme_mesh <= MDLtime_ext[-1]))
                    ext_mdl_oi = np.logical_and(MDLhgt_ext_mesh >= MEAShgt[0], MDLhgt_ext_mesh <= MAXhgt)
                    ###############################################################
                    # Set ranges based on the data instead of being hardcoded
                    #TEMPmin_ext,TEMPmax_ext,TEMPtick_ext = atmo_ylim((MEASt[ext_meas_oi],\
                    #                MDLt_ext[ext_mdl_oi]),plt_param.TEMPrnd_bse,plt_param.TICKint)    
                    #TEMPclabels_ext=np.arange(TEMPmin_ext,TEMPmax_ext,plt_param.TEMPcint)
                    #MIXRmin_ext,MIXRmax_ext,MIXRtick_ext = atmo_ylim((MEASmx[ext_meas_oi],\
                    #                MDLmx_ext[ext_mdl_oi]),plt_param.MIXRrnd_bse,plt_param.TICKint)    
                    #MIXRclabels_ext=np.arange(MIXRmin_ext,MIXRmax_ext,plt_param.MIXRcint)
                    ###############################################################
                    #Define Xtick Information (Push Out to Maximum Forecast Time)
                    Xmin_ext = MDLtime_ext[0];
                    Xmax_ext = (MDLtime_ext[0]+EXTlen*HRsec);
                    # do not use indice approach because we plan to push out fxcst 
                    ###############################################################
                    Xtick_ext = np.arange(MDLtime_ext[0],MDLtime_ext[0]+EXTlen*HRsec+1e-5,plt_param.Xdel_ext)
                    Xtick_ext_lab = [time.strftime('%H', time.gmtime(time_ref)) for time_ref in Xtick_ext]
                    Time_range_dates_ext = [basedate+timedelta(seconds=(time_ref)) for time_ref in Xtick_ext]
                    Date_start_label_ext=Time_range_dates_ext[0].strftime('%Y-%m-%d')
                    Date_end_label_ext=Time_range_dates_ext[-1].strftime('%Y-%m-%d')
                    ###############################################################
                
                    """ Define Bounds for the Model-Specific Standard Forecast Period """
                    ###############################################################
                    # Adjust Model Variables #
                    MDLt = MDLt[:,0:FXlen]; MDLmx = MDLmx[:,0:FXlen];
                    MDLtheta = MDLtheta[:,0:FXlen];
                    MDLMEASt = MDLMEASt[:,0:FXlen]; MDLMEASmx = MDLMEASmx[:,0:FXlen];
                    MDLMEAStheta = MDLMEAStheta[:,0:FXlen];
                    MDLtime = MDLtime[0:FXlen];
                    ###############################################################
                    """ Perform Model Averaging for Plotting and Statistical Analyses """
                    MDLtme_mesh, MDLhgt_mesh = np.meshgrid(MDLtime,MDLhgt); #define reference meshgrid
                    MDLMEAStme_mesh, MDLMEAShgt_mesh = np.meshgrid(MDLtime,MEAShgt)
                    ###############################################################
                    if (version == 'ASSIST'):
                        MEASt_mu = atmo_meas2grid(MEASt_nan,(MEAShgt_mesh,MEAStme_mesh),(MEAShgt,MDLtime),AVGmins=AVGper)
                        MEASmx_mu = atmo_meas2grid(MEASmx_nan,(MEAShgt_mesh,MEAStme_mesh),(MEAShgt,MDLtime),AVGmins=AVGper)
                        MEAStheta_mu = atmo_meas2grid(MEAStheta_nan,(MEAShgt_mesh,MEAStme_mesh),(MEAShgt,MDLtime),AVGmins=AVGper)
                    else:
                        MEASt_mu = atmo_meas2grid(MEASt,(MEAShgt_mesh,MEAStme_mesh),(MEAShgt,MDLtime),AVGmins=AVGper)
                        MEASmx_mu = atmo_meas2grid(MEASmx,(MEAShgt_mesh,MEAStme_mesh),(MEAShgt,MDLtime),AVGmins=AVGper)
                        MEAStheta_mu = atmo_meas2grid(MEAStheta,(MEAShgt_mesh,MEAStme_mesh),(MEAShgt,MDLtime),AVGmins=AVGper)
                    ###############################################################
                    meas_oi = np.logical_and(MEAShgt_mesh <= MAXhgt,np.logical_and(MEAStme_mesh >= MDLtime[0], MEAStme_mesh <= MDLtime[FXlen-1]))
                    mdl_oi = np.logical_and(MDLhgt_mesh >= MEAShgt[0], MDLhgt_mesh <= MAXhgt)
                    ###############################################################
                    # Set ranges based on the data instead of being hardcoded
                    #TEMPmin,TEMPmax,TEMPtick = atmo_ylim((MEASt[meas_oi],MDLt[mdl_oi]),\
                    #                        plt_param.TEMPrnd_bse,plt_param.TICKint)    
                    #TEMPclabels=np.arange(TEMPmin,TEMPmax,plt_param.TEMPcint)
                    #MIXRmin,MIXRmax,MIXRtick = atmo_ylim((MEASmx[meas_oi],MDLmx[mdl_oi]),\
                    #                        plt_param.MIXRrnd_bse,plt_param.TICKint)    
                    #MIXRclabels=np.arange(MIXRmin,MIXRmax,plt_param.MIXRcint)
                    ###############################################################
                    #Define Xtick Information (Push Out to Maximum Forecast Time)
                    Xmin = MDLtime[0];
                    Xmax = MDLtime[0]+(FXlen-1)*HRsec;
                    # do not use indice approach because we plan to push out fxcst 
                    ###############################################################
                    Xtick = np.arange(MDLtime[0],MDLtime[0]+(FXlen-1)*HRsec+1e-5,plt_param.Xdel)
                    Xtick_lab = [time.strftime('%H', time.gmtime(time_ref)) for time_ref in Xtick]
                    Time_range_dates = [basedate+timedelta(seconds=(time_ref)) for time_ref in Xtick]
                    Date_start_label=Time_range_dates[0].strftime('%Y-%m-%d')
                    Date_end_label=Time_range_dates[-1].strftime('%Y-%m-%d')
                    
                else:
                    """ Define Bounds for the Model-Specific Standard Forecast Period """
                    ###############################################################
                    # Know that an Extended Forecast is Non-Existent. Therefore,
                    # Simply Define meas_oi based on the Desired Forecast Time of 
                    # 86400
                    ###############################################################
                    """ Perform Model Averaging for Plotting and Statistical Analyses """
                    MDLtme_mesh, MDLhgt_mesh = np.meshgrid(MDLtime,MDLhgt); #define reference meshgrid
                    MDLMEAStme_mesh, MDLMEAShgt_mesh = np.meshgrid(MDLtime,MEAShgt)
                    ###############################################################
                    meas_oi = np.logical_and(MEAShgt_mesh <= MAXhgt,np.logical_and(MEAStme_mesh >= MDLtime[0], MEAStme_mesh <= MDLtime[0]+(FXlen-1)*HRsec))
                    mdl_oi = np.logical_and(MDLhgt_mesh >= MEAShgt[0], MDLhgt_mesh <= MAXhgt)
                    ###############################################################
                    if (version == 'ASSIST'):
                        MEASt_mu = atmo_meas2grid(MEASt_nan,(MEAShgt_mesh,MEAStme_mesh),(MEAShgt,MDLtime),AVGmins=AVGper)
                        MEASmx_mu = atmo_meas2grid(MEASmx_nan,(MEAShgt_mesh,MEAStme_mesh),(MEAShgt,MDLtime),AVGmins=AVGper)
                        MEAStheta_mu = atmo_meas2grid(MEAStheta_nan,(MEAShgt_mesh,MEAStme_mesh),(MEAShgt,MDLtime),AVGmins=AVGper)
                    else:
                        MEASt_mu = atmo_meas2grid(MEASt,(MEAShgt_mesh,MEAStme_mesh),(MEAShgt,MDLtime),AVGmins=AVGper)
                        MEASmx_mu = atmo_meas2grid(MEASmx,(MEAShgt_mesh,MEAStme_mesh),(MEAShgt,MDLtime),AVGmins=AVGper)
                        MEAStheta_mu = atmo_meas2grid(MEAStheta,(MEAShgt_mesh,MEAStme_mesh),(MEAShgt,MDLtime),AVGmins=AVGper)
                    ###############################################################
                    # Set ranges based on the data instead of being hardcoded
                    #TEMPmin,TEMPmax,TEMPtick = atmo_ylim((MEASt[meas_oi],MDLt[mdl_oi]),\
                    #                        plt_param.TEMPrnd_bse,plt_param.TICKint)    
                    #TEMPclabels=np.arange(TEMPmin,TEMPmax,plt_param.TEMPcint)
                    #MIXRmin,MIXRmax,MIXRtick = atmo_ylim((MEASmx[meas_oi],MDLmx[mdl_oi]),\
                    #                        plt_param.MIXRrnd_bse,plt_param.TICKint)    
                    #MIXRclabels=np.arange(MIXRmin,MIXRmax,plt_param.MIXRcint)
                    ###############################################################
                    #Define Xtick Information (Push Out to Maximum Forecast Time)
                    Xmin = MDLtime[0];
                    Xmax = MDLtime[0]+(FXlen-1)*HRsec;
                    # do not use indice approach because we plan to push out fxcst 
                    ###############################################################
                    Xtick = np.arange(MDLtime[0],MDLtime[0]+(FXlen-1)*HRsec+1e-5,plt_param.Xdel)
                    Xtick_lab = [time.strftime('%H', time.gmtime(time_ref)) for time_ref in Xtick]
                    Time_range_dates = [basedate+timedelta(seconds=(time_ref)) for time_ref in Xtick]
                    Date_start_label=Time_range_dates[0].strftime('%Y-%m-%d')
                    Date_end_label=Time_range_dates[-1].strftime('%Y-%m-%d')
                    
    
                """ ###############################################################
                Temperature Plots (instruct to share x-axis)
                ############################################################### """
                temp_fig, (temp_ax,mdl_ax,bias_ax) = plt.subplots(3,1,sharex = True,figsize=(plt_param.FIGwdth,plt_param.FIGwdth),\
                         height_ratios=[2,2,1],\
                         facecolor = 'w',edgecolor ='w',clear='True')
                ###################################################################
                # Plot the Temperature Time Histories
                ##########################################
                #temp_fig.suptitle('TROPoe '+version+' TEMP',x = plt_param.DLmeas_x_title,\
                temp_fig.suptitle(' TEMP',x = plt_param.DLmeas_x_title,\
                                y = plt_param.DLmeas_y_title2,horizontalalignment='right',\
                        fontsize=plt_param.TITLEfs);
                ###################################################################
                if (version == 'ASSIST'):
                    temp = temp_ax.pcolormesh(MEAStme_mesh,MEAShgt_mesh,MEASt,vmin = TEMPmin,\
                                            vmax = TEMPmax,shading='nearest', cmap = cmap_temp, alpha=0.5) 
                    temp2 = temp_ax.pcolormesh(MEAStme_mesh,MEAShgt_mesh,MEASt_nan,vmin = TEMPmin,\
                                            vmax = TEMPmax,shading='nearest', cmap = cmap_temp) 
                else:
                    temp2 = temp_ax.pcolormesh(MEAStme_mesh,MEAShgt_mesh,MEASt,vmin = TEMPmin,\
                                            vmax = TEMPmax,shading='nearest', cmap = cmap_temp) 
                temp_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5);
                ###################################################################
                cbar = temp_fig.colorbar(temp2,ax=temp_ax,ticks = TEMPtick,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Temperature, $\mathrm{deg}^{\circ}$ C', fontsize=plt_param.LABELfs);
                ###################################################################
                # overlay cloud base height 
                temp_ax.plot(MEAStme_mesh[0],MEAScbh_nondef,'o',markeredgewidth=plt_param.PLTmewidth,\
                    markerfacecolor=plt_param.PLTmfcolor,markeredgecolor=plt_param.PLTmecolor,\
                    markersize=plt_param.PLTmsize,label=GMLcbhlabel)
                temp_ax.plot(MEAStme_mesh[0],MEAScbh_def,'o',markeredgewidth=plt_param.PLTmewidth,\
                    markerfacecolor=plt_param.PLTmfcolor_miss,markeredgecolor=plt_param.PLTmecolor,\
                    markersize=plt_param.PLTmsize,label=DEFcbhlabel)
                temp_ax.legend(loc=8,fontsize=plt_param.LEGfs_bottom,bbox_to_anchor = plt_param.LEGpad_extra_bottom,\
                             handletextpad=.2, markerscale = 1, ncol = 2,frameon = False);
                # overlay contours
                tc = temp_ax.contour(MEAStme_mesh,MEAShgt_mesh,MEASt,levels=TEMPclabels,colors=ccolor,\
                    linewidths=plt_param.PLTcwidth,linestyles='solid')
                cb = temp_ax.clabel(tc,colors=ccolor,inline=True,inline_spacing=-1,fontsize=plt_param.PLTclsize)
                [cl.set_rotation(0) for cl in cb]
                ###################################################################
                temp_ax.set_ylim(HGTmin,HGTmax); temp_ax.set_yticks(HGTtick);temp_ax.set_yticklabels(HGTtick_label);
                temp_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                # set figure y-label
                temp_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                temp_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                ###################################################################
                #insert data date as title
                T = temp_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                T.set_position(plt_param.Tloc);
                ###############################
    
                ###################################################################
                # Plot the Model Temperature Time Histories
                ##########################################
                temp_mdl = mdl_ax.pcolormesh(MDLtme_mesh,MDLhgt_mesh,MDLt,vmin = TEMPmin,\
                                            vmax = TEMPmax,shading='gouraud', cmap = cmap_temp)   
                mdl_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                cbar = temp_fig.colorbar(temp_mdl,ax=mdl_ax,ticks = TEMPtick,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Temperature, $\mathrm{deg}^{\circ}$ C', fontsize=plt_param.LABELfs);
                ###################################################################
                # overlay contours
                if (MDLt.shape[1] > 1):
                    tc = mdl_ax.contour(MDLtme_mesh,MDLhgt_mesh,MDLt,levels=TEMPclabels,colors=ccolor,\
                        linewidths=plt_param.PLTcwidth,linestyles='solid')
                    cb = mdl_ax.clabel(tc,colors=ccolor,inline=True,inline_spacing=-1,fontsize=plt_param.PLTclsize)
                    [cl.set_rotation(0) for cl in cb]
                ###################################################################
                mdl_ax.set_ylim(HGTmin,HGTmax); mdl_ax.set_yticks(HGTtick);mdl_ax.set_yticklabels(HGTtick_label);
                mdl_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                ###################################################################
                # set figure y-label
                mdl_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                mdl_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                ###################################################################
                # Denote that Model Data was Limited
                if len(MDLtime) == 1:
                    mdl_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([HGTmin,MAXhgt]),'Single Forecast Hour Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
                # plotting function will auto-expand the grid
                ###################################################################
                M = mdl_ax.set_title(MDLname + ' TEMP'+MDLstr,\
                                     {'fontsize': plt_param.TITLEfs},pad=1,\
                             horizontalalignment='right');  
                M.set_position(plt_param.SPTLt_loc);
                
                """ ###############################################################
                Define and Plot the Model Temperature Biases
                ############################################################### """
                TEMPmdl_bias = MDLMEASt - MEASt_mu;
                #######################################
                # Determine Appropriate Y-Lims for Bias
                #TEMPb_min, TEMPb_max,TEMPb_tick = atmo_sptlbias_ylim((TEMPmdl_bias),plt_param.TEMPsptl_bias_bse,plt_param.TICKint)
                ###################################################################
                # Plot the Model Temperature Bias
                ##########################################
                norm=colors.BoundaryNorm(boundaries=TEMPb_bounds,ncolors=256)
                #bias_temp = bias_ax.pcolormesh(MDLMEAStme_mesh,MDLMEAShgt_mesh,TEMPmdl_bias,vmin = TEMPb_min,\
                #                            vmax = TEMPb_max,shading='nearest',cmap = cmap_bias_temp)   
                bias_temp = bias_ax.pcolormesh(MDLMEAStme_mesh,MDLMEAShgt_mesh,TEMPmdl_bias,\
                                            shading='nearest',cmap = cmap_bias_temp, norm=norm)   
                bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                cbar = temp_fig.colorbar(bias_temp,ax=bias_ax,ticks = TEMPb_tick,pad=plt_param.Cbar_pad,\
                    aspect = plt_param.CBaspect, drawedges=True)
                cbar.ax.minorticks_off()
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Temp, $\mathrm{deg}^{\circ}$ C', fontsize=plt_param.LABELfs);
                ###################################################################
                bias_ax.set_ylim(HGTmin,HGTmax); bias_ax.set_yticks(HGTtick);bias_ax.set_yticklabels(HGTtick_label);
                bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                ###################################################################
                # set figure y-label
                bias_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                ###################################################################
                bias_ax.set_xticks(Xtick); bias_ax.set_xticklabels(Xtick_lab,fontsize=plt_param.TICKfs);
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=UserWarning)
                    bias_ax.set_xlim(Xmin,Xmax)
                ################################################################### 
                bias_ax.set_xlabel('Hour, UTC',{'fontsize': plt_param.LABELfs});
                bias_ax.xaxis.set_label_coords(plt_param.Xlab_pad_tri[0],plt_param.Xlab_pad_tri[1]);
                ###################################################################
                B = bias_ax.set_title('TEMP Error (Mdl - Obs)',{'fontsize': plt_param.TITLEfs},pad=1,\
                             horizontalalignment='right');  
                B.set_position(plt_param.SPTLt_loc);
                ###################################################################
                temp_fig.subplots_adjust(left = 0.1,right = 0.96,bottom = 0.125, top = 0.95)
                ###################################################################
                bias_ax.annotate(Date_start_label,xy = plt_param.DATEpad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                bias_ax.annotate(Date_end_label,xy = plt_param.DATE_end_pad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                #######################################################################
                # Add some extra info
                bias_ax.annotate(MDLgridinfo,xy = plt_param.GRIDpad_left_tri,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'left',\
                                fontsize = plt_param.TITLEfs)
                bias_ax.annotate(ATTRstr,xy = plt_param.ATTRpad_ctr_dual,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'center',\
                                fontsize = plt_param.TITLEfs)
                #######################################################################
                # add a logo
                imax = temp_fig.add_axes(plt_param.LOGOaxcoords)
                imax.set_axis_off()
                imax.imshow(logo, aspect="equal")
                ###################################################################
                #Save the Developed Figure
                FILout = os.path.join(IMGpth,FILpfx+MDLname+'_TEMP_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                temp_fig.savefig(FILout,dpi=300); plt.close('all');
                                     
    
                """ ###############################################################
                Mixing Ratio Plots (instruct to share x-axis)
                ############################################################### """
                mx_fig, (mx_ax,mdl_ax,bias_ax) = plt.subplots(3,1,sharex = True,figsize=(plt_param.FIGwdth,plt_param.FIGwdth),\
                         height_ratios=[2,2,1],\
                         facecolor = 'w',edgecolor ='w',clear='True')
                ###################################################################
                # Plot the Mixing Ratio Time Histories
                ##########################################
                #mx_fig.suptitle('TROPoe '+version+' MIXR',x = plt_param.DLmeas_x_title,\
                mx_fig.suptitle(' MIXR',x = plt_param.DLmeas_x_title,\
                                y = plt_param.DLmeas_y_title2,horizontalalignment='right',\
                        fontsize=plt_param.TITLEfs);
                ###################################################################
                if (version == 'ASSIST'):
                    mx = mx_ax.pcolormesh(MEAStme_mesh,MEAShgt_mesh,MEASmx,vmin = MIXRmin,\
                                            vmax = MIXRmax,shading='nearest', cmap = cmap_mx, alpha=0.5) 
                    mx2 = mx_ax.pcolormesh(MEAStme_mesh,MEAShgt_mesh,MEASmx_nan,vmin = MIXRmin,\
                                            vmax = MIXRmax,shading='nearest', cmap = cmap_mx) 
                else:
                    mx2 = mx_ax.pcolormesh(MEAStme_mesh,MEAShgt_mesh,MEASmx,vmin = MIXRmin,\
                                            vmax = MIXRmax,shading='nearest', cmap = cmap_mx) 
                mx_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5);
                ###################################################################
                cbar = mx_fig.colorbar(mx2,ax=mx_ax,ticks = MIXRtick,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Mixing Ratio, g $\mathrm{kg}^{-1}$', fontsize=plt_param.LABELfs);
                ###################################################################
                # overlay cloud base height 
                mx_ax.plot(MEAStme_mesh[0],MEAScbh_nondef,'o',markeredgewidth=plt_param.PLTmewidth,\
                    markerfacecolor=plt_param.PLTmfcolor,markeredgecolor=plt_param.PLTmecolor,\
                    markersize=plt_param.PLTmsize,label=GMLcbhlabel)
                mx_ax.plot(MEAStme_mesh[0],MEAScbh_def,'o',markeredgewidth=plt_param.PLTmewidth,\
                    markerfacecolor=plt_param.PLTmfcolor_miss,markeredgecolor=plt_param.PLTmecolor,\
                    markersize=plt_param.PLTmsize,label=DEFcbhlabel)
                mx_ax.legend(loc=8,fontsize=plt_param.LEGfs_bottom,bbox_to_anchor = plt_param.LEGpad_extra_bottom,\
                             handletextpad=.2, markerscale = 1, ncol = 2,frameon = False);
                # overlay contours
                tc = mx_ax.contour(MEAStme_mesh,MEAShgt_mesh,MEASmx,levels=MIXRclabels,colors=ccolor,\
                    linewidths=plt_param.PLTcwidth,linestyles='solid')
                cb = mx_ax.clabel(tc,colors=ccolor,inline=True,inline_spacing=-1,fontsize=plt_param.PLTclsize)
                [cl.set_rotation(0) for cl in cb]
                ###################################################################
                mx_ax.set_ylim(HGTmin,HGTmax); mx_ax.set_yticks(HGTtick);mx_ax.set_yticklabels(HGTtick_label);
                mx_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                # set figure y-label
                mx_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                mx_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                ###################################################################
                #insert data date as title
                T = mx_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                T.set_position(plt_param.Tloc);
                ###############################
    
                ###################################################################
                # Plot the Model Mixing Ratio Time Histories
                ##########################################
                mx_mdl = mdl_ax.pcolormesh(MDLtme_mesh,MDLhgt_mesh,MDLmx,vmin = MIXRmin,\
                                            vmax = MIXRmax,shading='gouraud', cmap = cmap_mx)   
                mdl_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                cbar = mx_fig.colorbar(mx_mdl,ax=mdl_ax,ticks = MIXRtick,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Mixing Ratio, g $\mathrm{kg}^{-1}$', fontsize=plt_param.LABELfs);
                ###################################################################
                # overlay contours
                if (MDLmx.shape[1] > 1):
                    tc = mdl_ax.contour(MDLtme_mesh,MDLhgt_mesh,MDLmx,levels=MIXRclabels,colors=ccolor,\
                        linewidths=plt_param.PLTcwidth,linestyles='solid')
                    cb = mdl_ax.clabel(tc,colors=ccolor,inline=True,inline_spacing=-1,fontsize=plt_param.PLTclsize)
                    [cl.set_rotation(0) for cl in cb]
                ###################################################################
                mdl_ax.set_ylim(HGTmin,HGTmax); mdl_ax.set_yticks(HGTtick);mdl_ax.set_yticklabels(HGTtick_label);
                mdl_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                ###################################################################
                # set figure y-label
                mdl_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                mdl_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                ###################################################################
                # Denote that Model Data was Limited
                if len(MDLtime) == 1:
                    mdl_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([HGTmin,MAXhgt]),'Single Forecast Hour Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
                # plotting function will auto-expand the grid
                ###################################################################
                M = mdl_ax.set_title(MDLname + ' MIXR'+MDLstr,\
                                     {'fontsize': plt_param.TITLEfs},pad=1,\
                             horizontalalignment='right');  
                M.set_position(plt_param.SPTLt_loc);
                
                """ ###############################################################
                Define and Plot the Model Mixing Ratio Biases
                ############################################################### """
                MIXRmdl_bias = MDLMEASmx - MEASmx_mu;
                #######################################
                # Determine Appropriate Y-Lims for Bias
                #MIXRb_min, MIXRb_max,MIXRb_tick = atmo_sptlbias_ylim((MIXRmdl_bias),plt_param.MIXRbias_bse,plt_param.TICKint)
                ###################################################################
                # Plot the Model Mixing Ratio Bias
                ##########################################
                norm=colors.BoundaryNorm(boundaries=MIXRb_bounds,ncolors=256)
                #bias_mx = bias_ax.pcolormesh(MDLMEAStme_mesh,MDLMEAShgt_mesh,MIXRmdl_bias,vmin = MIXRb_min,\
                #                            vmax = MIXRb_max,shading='nearest',cmap = cmap_bias_mx)   
                bias_mx = bias_ax.pcolormesh(MDLMEAStme_mesh,MDLMEAShgt_mesh,MIXRmdl_bias,\
                                            shading='nearest',cmap = cmap_bias_mx, norm=norm)   
                bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                cbar = mx_fig.colorbar(bias_mx,ax=bias_ax,ticks = MIXRb_tick,pad=plt_param.Cbar_pad,\
                    aspect = plt_param.CBaspect,drawedges=True)
                cbar.ax.minorticks_off()
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Mix Ratio, g $\mathrm{kg}^{-1}$', fontsize=plt_param.LABELfs);
                ###################################################################
                bias_ax.set_ylim(HGTmin,HGTmax); bias_ax.set_yticks(HGTtick);bias_ax.set_yticklabels(HGTtick_label);
                bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                ###################################################################
                # set figure y-label
                bias_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                ###################################################################
                bias_ax.set_xticks(Xtick); bias_ax.set_xticklabels(Xtick_lab,fontsize=plt_param.TICKfs);
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=UserWarning)
                    bias_ax.set_xlim(Xmin,Xmax)
                ################################################################### 
                bias_ax.set_xlabel('Hour, UTC',{'fontsize': plt_param.LABELfs});
                bias_ax.xaxis.set_label_coords(plt_param.Xlab_pad_tri[0],plt_param.Xlab_pad_tri[1]);
                ###################################################################
                B = bias_ax.set_title('MIXR Error (Mdl - Obs)',{'fontsize': plt_param.TITLEfs},pad=1,\
                             horizontalalignment='right');  
                B.set_position(plt_param.SPTLt_loc);
                ###################################################################
                mx_fig.subplots_adjust(left = 0.1,right = 0.96,bottom = 0.125, top = 0.95)
                ###################################################################
                bias_ax.annotate(Date_start_label,xy = plt_param.DATEpad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                bias_ax.annotate(Date_end_label,xy = plt_param.DATE_end_pad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                #######################################################################
                # Add some extra info
                bias_ax.annotate(MDLgridinfo,xy = plt_param.GRIDpad_left_tri,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'left',\
                                fontsize = plt_param.TITLEfs)
                bias_ax.annotate(ATTRstr,xy = plt_param.ATTRpad_ctr_dual,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'center',\
                                fontsize = plt_param.TITLEfs)
                #######################################################################
                # add a logo
                imax = mx_fig.add_axes(plt_param.LOGOaxcoords)
                imax.set_axis_off()
                imax.imshow(logo, aspect="equal")
                ###################################################################
                #Save the Developed Figure
                FILout = os.path.join(IMGpth,FILpfx+MDLname+'_MIXR_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                mx_fig.savefig(FILout,dpi=300); plt.close('all');
    
                """ ###############################################################
                Potential Temperature Plots (instruct to share x-axis)
                ############################################################### """
                theta_fig, (theta_ax,mdl_ax,bias_ax) = plt.subplots(3,1,sharex = True,figsize=(plt_param.FIGwdth,plt_param.FIGwdth),\
                         height_ratios=[2,2,1],\
                         facecolor = 'w',edgecolor ='w',clear='True')
                ###################################################################
                # Plot the Potential Temperature Time Histories
                ##########################################
                #theta_fig.suptitle('TROPoe '+version+' THETA',x = plt_param.DLmeas_x_title,\
                theta_fig.suptitle(' THETA',x = plt_param.DLmeas_x_title,\
                                y = plt_param.DLmeas_y_title2,horizontalalignment='right',\
                        fontsize=plt_param.TITLEfs);
                ###################################################################
                if (version == 'ASSIST'):
                    theta = theta_ax.pcolormesh(MEAStme_mesh,MEAShgt_mesh,MEAStheta,vmin = THETAmin,\
                                            vmax = THETAmax,shading='nearest', cmap = cmap_theta, alpha=0.5) 
                    theta2 = theta_ax.pcolormesh(MEAStme_mesh,MEAShgt_mesh,MEAStheta_nan,vmin = THETAmin,\
                                            vmax = THETAmax,shading='nearest', cmap = cmap_theta) 
                else:
                    theta2 = theta_ax.pcolormesh(MEAStme_mesh,MEAShgt_mesh,MEAStheta,vmin = THETAmin,\
                                            vmax = THETAmax,shading='nearest', cmap = cmap_theta) 
                theta_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5);
                ###################################################################
                cbar = theta_fig.colorbar(theta2,ax=theta_ax,ticks = THETAtick,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Potential Temperature, K', fontsize=plt_param.LABELfs);
                ###################################################################
                # overlay cloud base height 
                theta_ax.plot(MEAStme_mesh[0],MEAScbh_nondef,'o',markeredgewidth=plt_param.PLTmewidth,\
                    markerfacecolor=plt_param.PLTmfcolor,markeredgecolor=plt_param.PLTmecolor,\
                    markersize=plt_param.PLTmsize,label=GMLcbhlabel)
                theta_ax.plot(MEAStme_mesh[0],MEAScbh_def,'o',markeredgewidth=plt_param.PLTmewidth,\
                    markerfacecolor=plt_param.PLTmfcolor_miss,markeredgecolor=plt_param.PLTmecolor,\
                    markersize=plt_param.PLTmsize,label=DEFcbhlabel)
                theta_ax.legend(loc=8,fontsize=plt_param.LEGfs_bottom,bbox_to_anchor = plt_param.LEGpad_extra_bottom,\
                             handletextpad=.2, markerscale = 1, ncol = 2,frameon = False);
                # overlay contours
                tc = theta_ax.contour(MEAStme_mesh,MEAShgt_mesh,MEAStheta,levels=THETAclabels,colors=ccolor,\
                    linewidths=plt_param.PLTcwidth,linestyles='solid')
                cb = theta_ax.clabel(tc,colors=ccolor,inline=True,inline_spacing=-1,fontsize=plt_param.PLTclsize)
                [cl.set_rotation(0) for cl in cb]
                ###################################################################
                theta_ax.set_ylim(HGTmin,HGTmax); theta_ax.set_yticks(HGTtick);theta_ax.set_yticklabels(HGTtick_label);
                theta_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                # set figure y-label
                theta_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                theta_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                ###################################################################
                #insert data date as title
                T = theta_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                T.set_position(plt_param.Tloc);
                ###############################
    
                ###################################################################
                # Plot the Model Potential Temperature Time Histories
                ##########################################
                theta_mdl = mdl_ax.pcolormesh(MDLtme_mesh,MDLhgt_mesh,MDLtheta,vmin = THETAmin,\
                                            vmax = THETAmax,shading='gouraud', cmap = cmap_theta)   
                mdl_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                cbar = theta_fig.colorbar(theta_mdl,ax=mdl_ax,ticks = THETAtick,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Potential Temperature, K', fontsize=plt_param.LABELfs);
                ###################################################################
                # overlay contours
                if (MDLtheta.shape[1] > 1):
                    tc = mdl_ax.contour(MDLtme_mesh,MDLhgt_mesh,MDLtheta,levels=THETAclabels,colors=ccolor,\
                        linewidths=plt_param.PLTcwidth,linestyles='solid')
                    cb = mdl_ax.clabel(tc,colors=ccolor,inline=True,inline_spacing=-1,fontsize=plt_param.PLTclsize)
                    [cl.set_rotation(0) for cl in cb]
                ###################################################################
                mdl_ax.set_ylim(HGTmin,HGTmax); mdl_ax.set_yticks(HGTtick);mdl_ax.set_yticklabels(HGTtick_label);
                mdl_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                ###################################################################
                # set figure y-label
                mdl_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                mdl_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                ###################################################################
                # Denote that Model Data was Limited
                if len(MDLtime) == 1:
                    mdl_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([HGTmin,MAXhgt]),'Single Forecast Hour Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
                # plotting function will auto-expand the grid
                ###################################################################
                M = mdl_ax.set_title(MDLname + ' THETA'+MDLstr,\
                                     {'fontsize': plt_param.TITLEfs},pad=1,\
                             horizontalalignment='right');  
                M.set_position(plt_param.SPTLt_loc);
                
                """ ###############################################################
                Define and Plot the Model Potential Temperature Biases
                ############################################################### """
                THETAmdl_bias = MDLMEAStheta - MEAStheta_mu;
                #######################################
                # Determine Appropriate Y-Lims for Bias
                #THETAb_min, THETAb_max,THETAb_tick = atmo_sptlbias_ylim((THETAmdl_bias),plt_param.THETAbias_bse,plt_param.TICKint)
                ###################################################################
                # Plot the Model Potential Temperature Bias
                ##########################################
                norm=colors.BoundaryNorm(boundaries=THETAb_bounds,ncolors=256)
                #bias_theta = bias_ax.pcolormesh(MDLMEAStme_mesh,MDLMEAShgt_mesh,THETAmdl_bias,vmin = THETAb_min,\
                #                            vmax = THETAb_max,shading='nearest',cmap = cmap_bias_theta)   
                bias_theta = bias_ax.pcolormesh(MDLMEAStme_mesh,MDLMEAShgt_mesh,THETAmdl_bias,\
                                            shading='nearest',cmap = cmap_bias_theta, norm=norm)   
                bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                cbar = theta_fig.colorbar(bias_theta,ax=bias_ax,ticks = THETAb_tick,pad=plt_param.Cbar_pad,\
                    aspect = plt_param.CBaspect,drawedges=True)
                cbar.ax.minorticks_off()
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Potential Temp, K', fontsize=plt_param.LABELfs);
                ###################################################################
                bias_ax.set_ylim(HGTmin,HGTmax); bias_ax.set_yticks(HGTtick);bias_ax.set_yticklabels(HGTtick_label);
                bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                ###################################################################
                # set figure y-label
                bias_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                ###################################################################
                bias_ax.set_xticks(Xtick); bias_ax.set_xticklabels(Xtick_lab,fontsize=plt_param.TICKfs);
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=UserWarning)
                    bias_ax.set_xlim(Xmin,Xmax)
                ################################################################### 
                bias_ax.set_xlabel('Hour, UTC',{'fontsize': plt_param.LABELfs});
                bias_ax.xaxis.set_label_coords(plt_param.Xlab_pad_tri[0],plt_param.Xlab_pad_tri[1]);
                ###################################################################
                B = bias_ax.set_title('THETA Error (Mdl - Obs)',{'fontsize': plt_param.TITLEfs},pad=1,\
                             horizontalalignment='right');  
                B.set_position(plt_param.SPTLt_loc);
                ###################################################################
                theta_fig.subplots_adjust(left = 0.1,right = 0.96,bottom = 0.125, top = 0.95)
                ###################################################################
                bias_ax.annotate(Date_start_label,xy = plt_param.DATEpad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                bias_ax.annotate(Date_end_label,xy = plt_param.DATE_end_pad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                #######################################################################
                # Add some extra info
                bias_ax.annotate(MDLgridinfo,xy = plt_param.GRIDpad_left_tri,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'left',\
                                fontsize = plt_param.TITLEfs)
                bias_ax.annotate(ATTRstr,xy = plt_param.ATTRpad_ctr_dual,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'center',\
                                fontsize = plt_param.TITLEfs)
                #######################################################################
                # add a logo
                imax = theta_fig.add_axes(plt_param.LOGOaxcoords)
                imax.set_axis_off()
                imax.imshow(logo, aspect="equal")
                ###################################################################
                #Save the Developed Figure
                FILout = os.path.join(IMGpth,FILpfx+MDLname+'_THETA_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                theta_fig.savefig(FILout,dpi=300); plt.close('all');
    
                """ Plot Extended Forecast Period if Applicable """
                ###################################################################
                if 'MDLtime_ext' in locals():
                    
                    """ Extended Temperature Plots """
                    temp_fig, (temp_ax,mdl_ax,bias_ax) = plt.subplots(3,1,sharex = True,figsize=(plt_param.FIGwdth,plt_param.FIGwdth),\
                             height_ratios=[2,2,1],\
                             facecolor = 'w',edgecolor ='w',clear='True')
                    ###############################################################
                    # Plot the Temperature Time Histories
                    ##########################################
                    #temp_fig.suptitle('TROPoe '+version+' TEMP',x = plt_param.DLmeas_x_title,\
                    temp_fig.suptitle(' TEMP',x = plt_param.DLmeas_x_title,\
                                    y = plt_param.DLmeas_y_title2,horizontalalignment='right',\
                            fontsize=plt_param.TITLEfs);
                    ###############################################################
                    if (version == 'ASSIST'):
                        temp = temp_ax.pcolormesh(MEAStme_mesh,MEAShgt_mesh,MEASt,vmin = TEMPmin_ext,\
                                            vmax = TEMPmax_ext,shading='nearest',cmap = cmap_temp, alpha=0.5) 
                        temp2 = temp_ax.pcolormesh(MEAStme_mesh,MEAShgt_mesh,MEASt_nan,vmin = TEMPmin_ext,\
                                            vmax = TEMPmax_ext,shading='nearest',cmap = cmap_temp) 
                    else:
                        temp2 = temp_ax.pcolormesh(MEAStme_mesh,MEAShgt_mesh,MEASt,vmin = TEMPmin_ext,\
                                            vmax = TEMPmax_ext,shading='nearest',cmap = cmap_temp) 
                    temp_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5);
                    ###############################################################
                    cbar = temp_fig.colorbar(temp2,ax=temp_ax,ticks = TEMPtick_ext,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                    cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                    cbar.set_label('Temperature, $\mathrm{deg}^{\circ}$ C', fontsize=plt_param.LABELfs);
                    ###############################################################
                    # overlay cloud base height 
                    temp_ax.plot(MEAStme_mesh[0],MEAScbh_nondef,'o',markeredgewidth=plt_param.PLTmewidth,\
                        markerfacecolor=plt_param.PLTmfcolor,markeredgecolor=plt_param.PLTmecolor,\
                        markersize=plt_param.PLTmsize,label=GMLcbhlabel)
                    temp_ax.plot(MEAStme_mesh[0],MEAScbh_def,'o',markeredgewidth=plt_param.PLTmewidth,\
                        markerfacecolor=plt_param.PLTmfcolor_miss,markeredgecolor=plt_param.PLTmecolor,\
                        markersize=plt_param.PLTmsize,label=DEFcbhlabel)
                    temp_ax.legend(loc=8,fontsize=plt_param.LEGfs_bottom,bbox_to_anchor = plt_param.LEGpad_extra_bottom,\
                                 handletextpad=.2, markerscale = 1, ncol = 2,frameon = False);
                    # overlay contours
                    tc = temp_ax.contour(MEAStme_mesh,MEAShgt_mesh,MEASt,levels=TEMPclabels_ext,colors=ccolor,\
                        linewidths=plt_param.PLTcwidth,linestyles='solid')
                    temp_ax.clabel(tc,colors=ccolor,inline=True,inline_spacing=-1,fontsize=plt_param.PLTclsize)
                    ###############################################################
                    temp_ax.set_ylim(HGTmin,HGTmax); temp_ax.set_yticks(HGTtick);temp_ax.set_yticklabels(HGTtick_label);
                    temp_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                    ###################################################################
                    # set figure y-label
                    temp_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                    temp_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                    ###############################################################
                    #insert data date as title
                    T = temp_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                    T.set_position(plt_param.Tloc);
                    ###############################
        
                    ###############################################################
                    # Plot the Model Temperature Time Histories
                    ##########################################
                    temp_mdl = mdl_ax.pcolormesh(MDLtme_ext_mesh,MDLhgt_ext_mesh,MDLt_ext,\
                                               vmin = TEMPmin_ext,vmax = TEMPmax_ext,shading='gouraud',cmap = cmap_temp)   
                    mdl_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                    ###############################################################
                    cbar = temp_fig.colorbar(temp_mdl,ax=mdl_ax,ticks = TEMPtick_ext,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                    cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                    cbar.set_label('Temperature, $\mathrm{deg}^{\circ}$ C', fontsize=plt_param.LABELfs);
                    ###############################################################
                    # overlay contours
                    if (MDLt_ext.shape[1] > 1):
                        tc = mdl_ax.contour(MDLtme_ext_mesh,MDLhgt_ext_mesh,MDLt_ext,levels=TEMPclabels,colors=ccolor,\
                            linewidths=plt_param.PLTcwidth,linestyles='solid')
                        mdl_ax.clabel(tc,colors=ccolor,inline=True,inline_spacing=-1,fontsize=plt_param.PLTclsize)
                    ###############################################################
                    mdl_ax.set_ylim(HGTmin,HGTmax); mdl_ax.set_yticks(HGTtick);mdl_ax.set_yticklabels(HGTtick_label);
                    mdl_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                    ###############################################################
                    # set figure y-label
                    mdl_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                    mdl_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                    ###############################################################
                    M = mdl_ax.set_title(MDLname + ' TEMP'+MDLstr,\
                            {'fontsize': plt_param.TITLEfs},pad=1,horizontalalignment='right');  
                    M.set_position(plt_param.SPTLt_loc);
                    
                    """ ###########################################################
                    Define and Plot the Model Temperature Biases
                    ########################################################### """
                    TEMPmdl_bias = MDLMEASt_ext - MEASt_ext_mu;
                    #######################################
                    # Determine Appropriate Y-Lims for Bias
                    #TEMPb_min, TEMPb_max,TEMPb_tick = atmo_sptlbias_ylim((TEMPmdl_bias),plt_param.TEMPsptl_bias_bse,plt_param.TICKint)
                    ###############################################################
                    # Plot the Model Temperature Bias
                    ##########################################
                    norm=colors.BoundaryNorm(boundaries=TEMPb_bounds,ncolors=256)
                    #bias_temp = bias_ax.pcolormesh(MDLMEAStme_ext_mesh,MDLMEAShgt_ext_mesh,\
                    #     TEMPmdl_bias,vmin = TEMPb_min,vmax = TEMPb_max,shading='nearest',cmap = cmap_bias_temp)
                    bias_temp = bias_ax.pcolormesh(MDLMEAStme_ext_mesh,MDLMEAShgt_ext_mesh,TEMPmdl_bias,\
                                            shading='nearest',cmap = cmap_bias_temp, norm=norm)   
                    bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                    ###############################################################
                    cbar = temp_fig.colorbar(bias_temp,ax=bias_ax,ticks = TEMPb_tick,pad=plt_param.Cbar_pad,\
                        aspect = plt_param.CBaspect,drawedges=True)
                    cbar.ax.minorticks_off()
                    cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                    cbar.set_label('Temp, $\mathrm{deg}^{\circ}$ C', fontsize=plt_param.LABELfs);
                    ###############################################################
                    bias_ax.set_ylim(HGTmin,HGTmax); bias_ax.set_yticks(HGTtick);bias_ax.set_yticklabels(HGTtick_label);
                    bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                    ###################################################################
                    # set figure y-label
                    bias_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                    bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                    ###############################################################
                    bias_ax.set_xticks(Xtick_ext); bias_ax.set_xticklabels(Xtick_ext_lab,fontsize=plt_param.TICKfs);
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", category=UserWarning)
                        bias_ax.set_xlim(Xmin_ext,Xmax_ext)
                    ############################################################### 
                    bias_ax.set_xlabel('Hour, UTC',{'fontsize': plt_param.LABELfs});
                    bias_ax.xaxis.set_label_coords(plt_param.Xlab_pad_tri[0],plt_param.Xlab_pad_tri[1]);
                    ###############################################################
                    B = bias_ax.set_title('TEMP Error (Mdl - Obs)',{'fontsize': plt_param.TITLEfs},pad=1,\
                                 horizontalalignment='right');  
                    B.set_position(plt_param.SPTLt_loc);
                    ###############################################################
                    temp_fig.subplots_adjust(left = 0.1,right = 0.96,bottom = 0.125, top = 0.95)
                    ###################################################################
                    bias_ax.annotate(Date_start_label_ext,xy = plt_param.DATEpad_sptl,xycoords = 'figure fraction',\
                                    verticalalignment = 'bottom',horizontalalignment = 'right',\
                                    fontsize = plt_param.TITLEfs)
                    bias_ax.annotate(Date_end_label_ext,xy = plt_param.DATE_end_pad_sptl,xycoords = 'figure fraction',\
                                    verticalalignment = 'bottom',horizontalalignment = 'right',\
                                    fontsize = plt_param.TITLEfs)
                    #######################################################################
                    # Add some extra info
                    bias_ax.annotate(MDLgridinfo,xy = plt_param.GRIDpad_left_tri,xycoords = 'figure fraction',\
                                    verticalalignment = 'bottom',horizontalalignment = 'left',\
                                    fontsize = plt_param.TITLEfs)
                    bias_ax.annotate(ATTRstr,xy = plt_param.ATTRpad_ctr_dual,xycoords = 'figure fraction',\
                                    verticalalignment = 'bottom',horizontalalignment = 'center',\
                                    fontsize = plt_param.TITLEfs)
                    #######################################################################
                    # add a logo
                    imax = temp_fig.add_axes(plt_param.LOGOaxcoords)
                    imax.set_axis_off()
                    imax.imshow(logo, aspect="equal")
                    ###############################################################          
                    #Save the Developed Figure
                    FILout = os.path.join(IMGpth_ext,FILpfx+MDLname+'_TEMP_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                    temp_fig.savefig(FILout,dpi=300); plt.close('all');
    
                    """ Extended Mixing Ratio Plots """
                    mx_fig, (mx_ax,mdl_ax,bias_ax) = plt.subplots(3,1,sharex = True,figsize=(plt_param.FIGwdth,plt_param.FIGwdth),\
                             height_ratios=[2,2,1],\
                             facecolor = 'w',edgecolor ='w',clear='True')
                    ###############################################################
                    # Plot the Mixing Ratio Time Histories
                    ##########################################
                    #mx_fig.suptitle('TROPoe '+version+' MIXR',x = plt_param.DLmeas_x_title,\
                    mx_fig.suptitle(' MIXR',x = plt_param.DLmeas_x_title,\
                                    y = plt_param.DLmeas_y_title2,horizontalalignment='right',\
                            fontsize=plt_param.TITLEfs);
                    ###############################################################
                    if (version == 'ASSIST'):
                        mx = mx_ax.pcolormesh(MEAStme_mesh,MEAShgt_mesh,MEASmx,vmin = MIXRmin_ext,\
                                            vmax = MIXRmax_ext,shading='nearest',cmap = cmap_mx,alpha=0.5) 
                        mx2 = mx_ax.pcolormesh(MEAStme_mesh,MEAShgt_mesh,MEASmx_nan,vmin = MIXRmin_ext,\
                                            vmax = MIXRmax_ext,shading='nearest',cmap = cmap_mx) 
                    else:
                        mx2 = mx_ax.pcolormesh(MEAStme_mesh,MEAShgt_mesh,MEASmx,vmin = MIXRmin_ext,\
                                            vmax = MIXRmax_ext,shading='nearest',cmap = cmap_mx) 
                    mx_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5);
                    ###############################################################
                    cbar = mx_fig.colorbar(mx2,ax=mx_ax,ticks = MIXRtick_ext,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                    cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                    cbar.set_label('Mixing Ratio, g $\mathrm{kg}^{-1}$', fontsize=plt_param.LABELfs);
                    ###############################################################
                    # overlay cloud base height 
                    mx_ax.plot(MEAStme_mesh[0],MEAScbh_nondef,'o',markeredgewidth=plt_param.PLTmewidth,\
                        markerfacecolor=plt_param.PLTmfcolor,markeredgecolor=plt_param.PLTmecolor,\
                        markersize=plt_param.PLTmsize,label=GMLcbhlabel)
                    mx_ax.plot(MEAStme_mesh[0],MEAScbh_def,'o',markeredgewidth=plt_param.PLTmewidth,\
                        markerfacecolor=plt_param.PLTmfcolor_miss,markeredgecolor=plt_param.PLTmecolor,\
                        markersize=plt_param.PLTmsize,label=DEFcbhlabel)
                    mx_ax.legend(loc=8,fontsize=plt_param.LEGfs_bottom,bbox_to_anchor = plt_param.LEGpad_extra_bottom,\
                                 handletextpad=.2, markerscale = 1, ncol = 2,frameon = False);
                    # overlay contours
                    tc = mx_ax.contour(MEAStme_mesh,MEAShgt_mesh,MEASmx,levels=MIXRclabels_ext,colors=ccolor,\
                        linewidths=plt_param.PLTcwidth,linestyles='solid')
                    mx_ax.clabel(tc,colors=ccolor,inline=True,inline_spacing=-1,fontsize=plt_param.PLTclsize)
                    ###############################################################
                    mx_ax.set_ylim(HGTmin,HGTmax); mx_ax.set_yticks(HGTtick);mx_ax.set_yticklabels(HGTtick_label);
                    mx_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                    ###################################################################
                    # set figure y-label
                    mx_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                    mx_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                    ###############################################################
                    #insert data date as title
                    T = mx_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                    T.set_position(plt_param.Tloc);
                    ###############################
        
                    ###############################################################
                    # Plot the Model Mixing Ratio Time Histories
                    ##########################################
                    mx_mdl = mdl_ax.pcolormesh(MDLtme_ext_mesh,MDLhgt_ext_mesh,MDLmx_ext,\
                                               vmin = MIXRmin_ext,vmax = MIXRmax_ext,shading='gouraud',cmap = cmap_mx)   
                    mdl_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                    ###############################################################
                    cbar = mx_fig.colorbar(mx_mdl,ax=mdl_ax,ticks = MIXRtick_ext,pad=plt_param.Cbar_pad,\
                            aspect = plt_param.CBaspect2)
                    cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                    cbar.set_label('Mixing Ratio, g $\mathrm{kg}^{-1}$', fontsize=plt_param.LABELfs);
                    ###############################################################
                    # overlay contours
                    if (MDLmx_ext.shape[1] > 1):
                        tc = mdl_ax.contour(MDLtme_ext_mesh,MDLhgt_ext_mesh,MDLmx_ext,levels=MIXRclabels,colors=ccolor,\
                            linewidths=plt_param.PLTcwidth,linestyles='solid')
                        mdl_ax.clabel(tc,colors=ccolor,inline=True,inline_spacing=-1,fontsize=plt_param.PLTclsize)
                    ###############################################################
                    mdl_ax.set_ylim(HGTmin,HGTmax); mdl_ax.set_yticks(HGTtick);mdl_ax.set_yticklabels(HGTtick_label);
                    mdl_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                    ###############################################################
                    # set figure y-label
                    mdl_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                    mdl_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                    ###############################################################
                    M = mdl_ax.set_title(MDLname + ' MIXR'+MDLstr,\
                            {'fontsize': plt_param.TITLEfs},pad=1,horizontalalignment='right');  
                    M.set_position(plt_param.SPTLt_loc);
                    
                    """ ###########################################################
                    Define and Plot the Model Mixing Ratio Biases
                    ########################################################### """
                    MIXRmdl_bias = MDLMEASmx_ext - MEASmx_ext_mu;
                    #######################################
                    # Determine Appropriate Y-Lims for Bias
                    #MIXRb_min, MIXRb_max,MIXRb_tick = atmo_sptlbias_ylim((MIXRmdl_bias),plt_param.MIXRbias_bse,plt_param.TICKint)
                    ###############################################################
                    # Plot the Model Mixing Ratio Bias
                    ##########################################
                    norm=colors.BoundaryNorm(boundaries=MIXRb_bounds,ncolors=256)
                    #bias_mx = bias_ax.pcolormesh(MDLMEAStme_ext_mesh,MDLMEAShgt_ext_mesh,\
                    #     MIXRmdl_bias,vmin = MIXRb_min,vmax = MIXRb_max,shading='nearest',cmap = cmap_bias_mx)
                    bias_mx = bias_ax.pcolormesh(MDLMEAStme_ext_mesh,MDLMEAShgt_ext_mesh,MIXRmdl_bias,\
                                                shading='nearest',cmap = cmap_bias_mx, norm=norm)   
                    bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                    ###############################################################
                    cbar = mx_fig.colorbar(bias_mx,ax=bias_ax,ticks = MIXRb_tick,pad=plt_param.Cbar_pad,\
                        aspect = plt_param.CBaspect,drawedges=True)
                    cbar.ax.minorticks_off()
                    cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                    cbar.set_label('Mix Ratio, g $\mathrm{kg}^{-1}$', fontsize=plt_param.LABELfs);
                    ###############################################################
                    bias_ax.set_ylim(HGTmin,HGTmax); bias_ax.set_yticks(HGTtick);bias_ax.set_yticklabels(HGTtick_label);
                    bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                    ###################################################################
                    # set figure y-label
                    bias_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                    bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                    ###############################################################
                    bias_ax.set_xticks(Xtick_ext); bias_ax.set_xticklabels(Xtick_ext_lab,fontsize=plt_param.TICKfs);
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", category=UserWarning)
                        bias_ax.set_xlim(Xmin_ext,Xmax_ext)
                    ############################################################### 
                    bias_ax.set_xlabel('Hour, UTC',{'fontsize': plt_param.LABELfs});
                    bias_ax.xaxis.set_label_coords(plt_param.Xlab_pad_tri[0],plt_param.Xlab_pad_tri[1]);
                    ###############################################################
                    B = bias_ax.set_title('MIXR Error (Mdl - Obs)',{'fontsize': plt_param.TITLEfs},pad=1,\
                                 horizontalalignment='right');  
                    B.set_position(plt_param.SPTLt_loc);
                    ###############################################################
                    mx_fig.subplots_adjust(left = 0.1,right = 0.96,bottom = 0.125, top = 0.95)
                    ###################################################################
                    bias_ax.annotate(Date_start_label_ext,xy = plt_param.DATEpad_sptl,xycoords = 'figure fraction',\
                                    verticalalignment = 'bottom',horizontalalignment = 'right',\
                                    fontsize = plt_param.TITLEfs)
                    bias_ax.annotate(Date_end_label_ext,xy = plt_param.DATE_end_pad_sptl,xycoords = 'figure fraction',\
                                    verticalalignment = 'bottom',horizontalalignment = 'right',\
                                    fontsize = plt_param.TITLEfs)
                    #######################################################################
                    # Add some extra info
                    bias_ax.annotate(MDLgridinfo,xy = plt_param.GRIDpad_left_tri,xycoords = 'figure fraction',\
                                    verticalalignment = 'bottom',horizontalalignment = 'left',\
                                    fontsize = plt_param.TITLEfs)
                    bias_ax.annotate(ATTRstr,xy = plt_param.ATTRpad_ctr_dual,xycoords = 'figure fraction',\
                                    verticalalignment = 'bottom',horizontalalignment = 'center',\
                                    fontsize = plt_param.TITLEfs)
                    #######################################################################
                    # add a logo
                    imax = mx_fig.add_axes(plt_param.LOGOaxcoords)
                    imax.set_axis_off()
                    imax.imshow(logo, aspect="equal")
                    ###############################################################          
                    #Save the Developed Figure
                    FILout = os.path.join(IMGpth_ext,FILpfx+MDLname+'_MIXR_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                    mx_fig.savefig(FILout,dpi=300); plt.close('all');
    
                    """ Extended Potential Temperature Plots """
                    theta_fig, (theta_ax,mdl_ax,bias_ax) = plt.subplots(3,1,sharex = True,figsize=(plt_param.FIGwdth,plt_param.FIGwdth),\
                             height_ratios=[2,2,1],\
                             facecolor = 'w',edgecolor ='w',clear='True')
                    ###############################################################
                    # Plot the Potential Temperature Time Histories
                    ##########################################
                    #theta_fig.suptitle('TROPoe '+version+' THETA',x = plt_param.DLmeas_x_title,\
                    theta_fig.suptitle(' THETA',x = plt_param.DLmeas_x_title,\
                                    y = plt_param.DLmeas_y_title2,horizontalalignment='right',\
                            fontsize=plt_param.TITLEfs);
                    ###############################################################
                    if (version == 'ASSIST'):
                        theta = theta_ax.pcolormesh(MEAStme_mesh,MEAShgt_mesh,MEAStheta,vmin = THETAmin_ext,\
                                            vmax = THETAmax_ext,shading='nearest',cmap = cmap_theta,alpha=0.5) 
                        theta2 = theta_ax.pcolormesh(MEAStme_mesh,MEAShgt_mesh,MEAStheta_nan,vmin = THETAmin_ext,\
                                            vmax = THETAmax_ext,shading='nearest',cmap = cmap_theta) 
                    else:
                        theta2 = theta_ax.pcolormesh(MEAStme_mesh,MEAShgt_mesh,MEAStheta,vmin = THETAmin_ext,\
                                            vmax = THETAmax_ext,shading='nearest',cmap = cmap_theta) 
                    theta_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5);
                    ###############################################################
                    cbar = theta_fig.colorbar(theta2,ax=theta_ax,ticks = THETAtick_ext,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                    cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                    cbar.set_label('Potential Temperature, K', fontsize=plt_param.LABELfs);
                    ###############################################################
                    # overlay cloud base height 
                    theta_ax.plot(MEAStme_mesh[0],MEAScbh_nondef,'o',markeredgewidth=plt_param.PLTmewidth,\
                        markerfacecolor=plt_param.PLTmfcolor,markeredgecolor=plt_param.PLTmecolor,\
                        markersize=plt_param.PLTmsize,label=GMLcbhlabel)
                    theta_ax.plot(MEAStme_mesh[0],MEAScbh_def,'o',markeredgewidth=plt_param.PLTmewidth,\
                        markerfacecolor=plt_param.PLTmfcolor_miss,markeredgecolor=plt_param.PLTmecolor,\
                        markersize=plt_param.PLTmsize,label='Cloud Base Height (default value, true value unknown)')
                    theta_ax.legend(loc=8,fontsize=plt_param.LEGfs_bottom,bbox_to_anchor = plt_param.LEGpad_extra_bottom,\
                                 handletextpad=.2, markerscale = 1, ncol = 2,frameon = False);
                    # overlay contours
                    tc = theta_ax.contour(MEAStme_mesh,MEAShgt_mesh,MEAStheta,levels=THETAclabels_ext,colors=ccolor,\
                        linewidths=plt_param.PLTcwidth,linestyles='solid')
                    theta_ax.clabel(tc,colors=ccolor,inline=True,inline_spacing=-1,fontsize=plt_param.PLTclsize)
                    ###############################################################
                    theta_ax.set_ylim(HGTmin,HGTmax); theta_ax.set_yticks(HGTtick);theta_ax.set_yticklabels(HGTtick_label);
                    theta_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                    ###################################################################
                    # set figure y-label
                    theta_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                    theta_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                    ###############################################################
                    #insert data date as title
                    T = theta_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                    T.set_position(plt_param.Tloc);
                    ###############################
        
                    ###############################################################
                    # Plot the Model Potential Temperature Time Histories
                    ##########################################
                    theta_mdl = mdl_ax.pcolormesh(MDLtme_ext_mesh,MDLhgt_ext_mesh,MDLtheta_ext,\
                                               vmin = THETAmin_ext,vmax = THETAmax_ext,shading='gouraud',cmap = cmap_theta)   
                    mdl_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                    ###############################################################
                    cbar = theta_fig.colorbar(theta_mdl,ax=mdl_ax,ticks = THETAtick_ext,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                    cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                    cbar.set_label('Potential Temperature, K', fontsize=plt_param.LABELfs);
                    ###############################################################
                    # overlay contours
                    if (MDLtheta_ext.shape[1] > 1):
                        tc = mdl_ax.contour(MDLtme_ext_mesh,MDLhgt_ext_mesh,MDLtheta_ext,levels=THETAclabels,colors=ccolor,\
                            linewidths=plt_param.PLTcwidth,linestyles='solid')
                        mdl_ax.clabel(tc,colors=ccolor,inline=True,inline_spacing=-1,fontsize=plt_param.PLTclsize)
                    ###############################################################
                    mdl_ax.set_ylim(HGTmin,HGTmax); mdl_ax.set_yticks(HGTtick);mdl_ax.set_yticklabels(HGTtick_label);
                    mdl_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                    ###############################################################
                    # set figure y-label
                    mdl_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                    mdl_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                    ###############################################################
                    M = mdl_ax.set_title(MDLname + ' THETA'+MDLstr,\
                            {'fontsize': plt_param.TITLEfs},pad=1,horizontalalignment='right');  
                    M.set_position(plt_param.SPTLt_loc);
                    
                    """ ###########################################################
                    Define and Plot the Model Potential Temperature Biases
                    ########################################################### """
                    THETAmdl_bias = MDLMEAStheta_ext - MEAStheta_ext_mu;
                    #######################################
                    # Determine Appropriate Y-Lims for Bias
                    #THETAb_min, THETAb_max,THETAb_tick = atmo_sptlbias_ylim((THETAmdl_bias),plt_param.THETAbias_bse,plt_param.TICKint)
                    ###############################################################
                    # Plot the Model Potential Temperature Bias
                    ##########################################
                    norm=colors.BoundaryNorm(boundaries=THETAb_bounds,ncolors=256)
                    #bias_theta = bias_ax.pcolormesh(MDLMEAStme_ext_mesh,MDLMEAShgt_ext_mesh,\
                    #     THETAmdl_bias,vmin = THETAb_min,vmax = THETAb_max,shading='nearest',cmap = cmap_bias_theta)
                    bias_theta = bias_ax.pcolormesh(MDLMEAStme_ext_mesh,MDLMEAShgt_ext_mesh,THETAmdl_bias,\
                                                shading='nearest',cmap = cmap_bias_theta, norm=norm)   
                    bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                    ###############################################################
                    cbar = theta_fig.colorbar(bias_theta,ax=bias_ax,ticks = THETAb_tick,pad=plt_param.Cbar_pad,\
                        aspect = plt_param.CBaspect,drawedges=True)
                    cbar.ax.minorticks_off()
                    cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                    cbar.set_label('Potential Temp, K', fontsize=plt_param.LABELfs);
                    ###############################################################
                    bias_ax.set_ylim(HGTmin,HGTmax); bias_ax.set_yticks(HGTtick);bias_ax.set_yticklabels(HGTtick_label);
                    bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                    ###################################################################
                    # set figure y-label
                    bias_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                    bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                    ###############################################################
                    bias_ax.set_xticks(Xtick_ext); bias_ax.set_xticklabels(Xtick_ext_lab,fontsize=plt_param.TICKfs);
                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore", category=UserWarning)
                        bias_ax.set_xlim(Xmin_ext,Xmax_ext)
                    ############################################################### 
                    bias_ax.set_xlabel('Hour, UTC',{'fontsize': plt_param.LABELfs});
                    bias_ax.xaxis.set_label_coords(plt_param.Xlab_pad_tri[0],plt_param.Xlab_pad_tri[1]);
                    ###############################################################
                    B = bias_ax.set_title('THETA Error (Mdl - Obs)',{'fontsize': plt_param.TITLEfs},pad=1,\
                                 horizontalalignment='right');  
                    B.set_position(plt_param.SPTLt_loc);
                    ###############################################################
                    theta_fig.subplots_adjust(left = 0.1,right = 0.96,bottom = 0.125, top = 0.95)
                    ###################################################################
                    bias_ax.annotate(Date_start_label_ext,xy = plt_param.DATEpad_sptl,xycoords = 'figure fraction',\
                                    verticalalignment = 'bottom',horizontalalignment = 'right',\
                                    fontsize = plt_param.TITLEfs)
                    bias_ax.annotate(Date_end_label_ext,xy = plt_param.DATE_end_pad_sptl,xycoords = 'figure fraction',\
                                    verticalalignment = 'bottom',horizontalalignment = 'right',\
                                    fontsize = plt_param.TITLEfs)
                    #######################################################################
                    # Add some extra info
                    bias_ax.annotate(MDLgridinfo,xy = plt_param.GRIDpad_left_tri,xycoords = 'figure fraction',\
                                    verticalalignment = 'bottom',horizontalalignment = 'left',\
                                    fontsize = plt_param.TITLEfs)
                    bias_ax.annotate(ATTRstr,xy = plt_param.ATTRpad_ctr_dual,xycoords = 'figure fraction',\
                                    verticalalignment = 'bottom',horizontalalignment = 'center',\
                                    fontsize = plt_param.TITLEfs)
                    #######################################################################
                    # add a logo
                    imax = theta_fig.add_axes(plt_param.LOGOaxcoords)
                    imax.set_axis_off()
                    imax.imshow(logo, aspect="equal")
                    ###############################################################          
                    #Save the Developed Figure
                    FILout = os.path.join(IMGpth_ext,FILpfx+MDLname+'_THETA_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                    theta_fig.savefig(FILout,dpi=300); plt.close('all');
    
                    del MDLtime_ext
    
            else:
                """ Initialization Time Not Available """
                ###################################################################
                MDLpseudo_time = np.arange(0,(FXlen-1)*HRsec+1e-5,step=HRsec)
                MDLtime = int(INIpos[ini])*HRsec + MDLpseudo_time
                MDLhgt = MEAShgt; #assumed model heights equal measurement heights
                ###################################################################
                #Define Xtick Information (Only Needs to be Defined Once)
                Xmin = MDLtime[0];
                Xmax = (MDLtime[0]+(FXlen-1)*HRsec);
                ###################################################################
                Xtick = np.arange(MDLtime[0],MDLtime[-1]+1e-5,plt_param.Xdel)
                Xtick_lab = [time.strftime('%H', time.gmtime(time_ref)) for time_ref in Xtick]
                Time_range_dates = [basedate+timedelta(seconds=(time_ref)) for time_ref in Xtick]
                Date_start_label=Time_range_dates[0].strftime('%Y-%m-%d')
                Date_end_label=Time_range_dates[-1].strftime('%Y-%m-%d')
    
                DLgrid = (MEAStime[0],MEAStime[-1],MEAShgt[0],MEAShgt[-1])     
                ###################################################################
                """ Define Relevant Parameters for Plotting """
                #MEAStme_base, MEAShgt_base = np.meshgrid(MEAStime,MEAShgt) #define reference meshgrid
                ###################################################################
                # Perform Correction to Create Equally Spaced Data Arrays #
                ###################################################################
                #MEAStme_mesh = MEAStme_base; MEAShgt_mesh = MEAShgt_base
                #MEAStme_mesh, MEAShgt_mesh, MEASt = atmo_sptl_correct(MEAStme_base,MEAShgt_base,MEASt);
                #_, _, MEAStd = atmo_sptl_correct(MEAStme_base,MEAShgt_base,MEAStd);
                #_, _, MEASrh = atmo_sptl_correct(MEAStme_base,MEAShgt_base,MEASrh);
                #_, _, MEASmx = atmo_sptl_correct(MEAStme_base,MEAShgt_base,MEASmx);
                #_, _, MEAScbh = atmo_sptl_correct(MEAStme_base,MEAShgt_base,MEAScbh);
                #MEAScbh_def = np.copy(MEAScbh)
                #MEAScbh_def[np.where(MEAScbh != DEFcbh)] = np.nan;
                #MEAScbh_nondef = np.copy(MEAScbh)
                #MEAScbh_nondef[np.where(MEAScbh==DEFcbh)] = np.nan;
                #if (version == 'ASSIST'):
                #    _, _, MEASt_nan = atmo_sptl_correct(MEAStme_base,MEAShgt_base,MEASt_nan);
                #    _, _, MEASmx_nan = atmo_sptl_correct(MEAStme_base,MEAShgt_base,MEASmx_nan);
    
                ###################################################################
                # Define Proper Indices to Determine Axes Parameters
                meas_oi = np.logical_and(MEAShgt_mesh <= MAXhgt,np.logical_and(MEAStme_mesh >= MDLtime[0], MEAStme_mesh <= MDLtime[-1]))
                ###################################################################
                # Set ranges based on the data instead of being hardcoded
                #TEMPmin,TEMPmax,TEMPtick = atmo_ylim((MEASt[meas_oi]),plt_param.TEMPrnd_bse,plt_param.TICKint)
                #TEMPclabels=np.arange(TEMPmin,TEMPmax,plt_param.TEMPcint)
                #MIXRmin,MIXRmax,MIXRtick = atmo_ylim((MEASt[meas_oi]),plt_param.MIXRrnd_bse,plt_param.TICKint) 
                #MIXRclabels=np.arange(MIXRmin,MIXRmax,plt_param.MIXRcint)
                ###################################################################
                # maximum and minimum measurement height 
                HGTmin,HGTmax,HGTtick = atmo_ylim((np.asarray([0,MAXhgt])),plt_param.SHGTrnd_bse,plt_param.TICKint);
                HGTtick_label = [str(ht/1000) for ht in HGTtick]
    
                """ ###############################################################
                Temperature (instruct to share x-axis)
                ############################################################### """
                temp_fig, (temp_ax,mdl_ax,bias_ax) = plt.subplots(3,1,sharex = True,figsize=(plt_param.FIGwdth,plt_param.FIGwdth),\
                         height_ratios=[2,2,1],\
                         facecolor = 'w',edgecolor ='w',clear='True')
                ###################################################################
                # Plot the Lidar Temperature Time Histories
                ##########################################
                #temp_fig.suptitle('TROPoe '+version+' TEMP',x = plt_param.DLmeas_x_title,\
                temp_fig.suptitle(' TEMP',x = plt_param.DLmeas_x_title,\
                                y = plt_param.DLmeas_y_title2,horizontalalignment='right',\
                        fontsize=plt_param.TITLEfs);
                ###################################################################
                #temp = temp_ax.pcolormesh(MEAStme_mesh,MEAShgt_mesh,MEASt,vmin = TEMPmin,\
                if (version == 'ASSIST'):
                    temp = temp_ax.pcolormesh(MEAStme_mesh,MEAShgt_mesh,MEASt,vmin = TEMPmin,\
                                            vmax = TEMPmax,shading='nearest', cmap = cmap_temp, alpha=0.5) 
                    temp2 = temp_ax.pcolormesh(MEAStme_mesh,MEAShgt_mesh,MEASt_nan,vmin = TEMPmin,\
                                            vmax = TEMPmax,shading='nearest', cmap = cmap_temp) 
                else:
                    temp2 = temp_ax.pcolormesh(MEAStme_mesh,MEAShgt_mesh,MEASt,vmin = TEMPmin,\
                                            vmax = TEMPmax,shading='nearest', cmap = cmap_temp) 
                temp_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5);
                ###################################################################
                cbar = temp_fig.colorbar(temp2,ax=temp_ax,ticks = TEMPtick,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Temperature, $\mathrm{deg}^{\circ}$ C', fontsize=plt_param.LABELfs);
                ###################################################################
                # overlay cloud base height 
                temp_ax.plot(MEAStme_mesh[0],MEAScbh_nondef,'o',markeredgewidth=plt_param.PLTmewidth,\
                    markerfacecolor=plt_param.PLTmfcolor,markeredgecolor=plt_param.PLTmecolor,\
                    markersize=plt_param.PLTmsize,label=GMLcbhlabel)
                temp_ax.plot(MEAStme_mesh[0],MEAScbh_def,'o',markeredgewidth=plt_param.PLTmewidth,\
                    markerfacecolor=plt_param.PLTmfcolor_miss,markeredgecolor=plt_param.PLTmecolor,\
                    markersize=plt_param.PLTmsize,label='Cloud Base Height (default value, true value unknown)')
                temp_ax.legend(loc=8,fontsize=plt_param.LEGfs_bottom,bbox_to_anchor = plt_param.LEGpad_extra_bottom,\
                             handletextpad=.2, markerscale = 1, ncol = 2,frameon = False);
                # overlay contours
                tc = temp_ax.contour(MEAStme_mesh,MEAShgt_mesh,MEASt,levels=TEMPclabels,colors=ccolor,\
                    linewidths=plt_param.PLTcwidth,linestyles='solid')
                temp_ax.clabel(tc,colors=ccolor,inline=True,inline_spacing=-1,fontsize=plt_param.PLTclsize)
                ###################################################################
                temp_ax.set_ylim(HGTmin,HGTmax); temp_ax.set_yticks(HGTtick);temp_ax.set_yticklabels(HGTtick_label);
                temp_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                # set figure y-label
                temp_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                temp_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                ###################################################################
                #insert data date as title
                T = temp_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                T.set_position(plt_param.Tloc);
                ###############################
    
                ###################################################################
                # Plot the Model Temperature Time Histories
                ##########################################
                temp_mdl = mdl_ax.imshow(MEASt*np.nan,extent=DLgrid,vmin=TEMPmin,vmax=TEMPmax,\
                                     origin='lower',aspect='auto', cmap=cmap_temp);   
                mdl_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                cbar = temp_fig.colorbar(temp_mdl,ax=mdl_ax,ticks = TEMPtick,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Temperature, $\mathrm{deg}^{\circ}$ C', fontsize=plt_param.LABELfs);
                ###################################################################
                mdl_ax.set_ylim(HGTmin,HGTmax); mdl_ax.set_yticks(HGTtick);mdl_ax.set_yticklabels(HGTtick_label);
                mdl_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                ###################################################################
                # set figure y-label
                mdl_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                mdl_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                ###################################################################
                # Denote that Model Data was Not Available
                mdl_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([HGTmin,MAXhgt]),MDLname+' Data Not Available',\
                         horizontalalignment='center',verticalalignment='bottom',color='red',\
                         fontweight='semibold',fontsize=plt_param.TITLEfs);
                ###################################################################
                M = mdl_ax.set_title(MDLname + ' TEMP'+MDLstr,\
                            {'fontsize': plt_param.TITLEfs},pad=1,horizontalalignment='right');  
                M.set_position(plt_param.SPTLt_loc);
                
                """ ###############################################################
                Define and Plot the Model Temperature Biases
                ############################################################### """
                # Plot the Model Temperature Bias
                ###################################################################
                #bias_temp = bias_ax.imshow(MEASt*np.nan,extent=DLgrid,vmin=TEMPb_min,vmax=TEMPb_max,\
                #                     origin='lower',aspect='auto', cmap=cmap_bias_temp);   
                norm=colors.BoundaryNorm(boundaries=TEMPb_bounds,ncolors=256)
                bias_temp = bias_ax.imshow(MEASt*np.nan,extent=DLgrid, \
                                     origin='lower',aspect='auto', cmap=cmap_bias_temp, norm=norm);   
                bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                cbar = temp_fig.colorbar(bias_temp,ax=bias_ax,ticks = TEMPb_tick,pad=plt_param.Cbar_pad,\
                    aspect = plt_param.CBaspect,drawedges=True)
                cbar.ax.minorticks_off()
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Temp, $\mathrm{deg}^{\circ}$ C', fontsize=plt_param.LABELfs);
                ###################################################################
                bias_ax.set_ylim(HGTmin,HGTmax); bias_ax.set_yticks(HGTtick);bias_ax.set_yticklabels(HGTtick_label);
                bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                ###################################################################
                # set figure y-label
                bias_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                ###################################################################
                bias_ax.set_xticks(Xtick); bias_ax.set_xticklabels(Xtick_lab,fontsize=plt_param.TICKfs);
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=UserWarning)
                    bias_ax.set_xlim(Xmin,Xmax)
                ###################################################################
                bias_ax.set_xlabel('Hour, UTC',{'fontsize': plt_param.LABELfs});
                bias_ax.xaxis.set_label_coords(plt_param.Xlab_pad_tri[0],plt_param.Xlab_pad_tri[1]);
                ###################################################################
                # Denote that Model Data was Not Available
                bias_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([HGTmin,MAXhgt]),MDLname+' Data Not Available',\
                         horizontalalignment='center',verticalalignment='bottom',color='red',\
                         fontweight='semibold',fontsize=plt_param.TITLEfs);
                ###################################################################
                B = bias_ax.set_title('TEMP Error (Mdl - Obs)',{'fontsize': plt_param.TITLEfs},pad=1,\
                             horizontalalignment='right');  
                B.set_position(plt_param.SPTLt_loc);
                ###################################################################
                temp_fig.subplots_adjust(left = 0.1,right = 0.96,bottom = 0.125, top = 0.95)
                ###################################################################
                bias_ax.annotate(Date_start_label,xy = plt_param.DATEpad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                bias_ax.annotate(Date_end_label,xy = plt_param.DATE_end_pad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                #######################################################################
                # Add some extra info
                bias_ax.annotate(MDLgridinfo,xy = plt_param.GRIDpad_left_tri,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'left',\
                                fontsize = plt_param.TITLEfs)
                bias_ax.annotate(ATTRstr,xy = plt_param.ATTRpad_ctr_dual,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'center',\
                                fontsize = plt_param.TITLEfs)
                #######################################################################
                # add a logo
                imax = temp_fig.add_axes(plt_param.LOGOaxcoords)
                imax.set_axis_off()
                imax.imshow(logo, aspect="equal")
                ###################################################################          
                #Save the Developed Figure
                FILout = os.path.join(IMGpth_model,FILpfx+MDLname+'_TEMP_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                temp_fig.savefig(FILout,dpi=300); plt.close('all');
                                     
                """ ###############################################################
                Mixing Ratio (instruct to share x-axis)
                ############################################################### """
                mx_fig, (mx_ax,mdl_ax,bias_ax) = plt.subplots(3,1,sharex = True,figsize=(plt_param.FIGwdth,plt_param.FIGwdth),\
                         height_ratios=[2,2,1],\
                         facecolor = 'w',edgecolor ='w',clear='True')
                ###################################################################
                # Plot the Lidar Mixing Ratio Time Histories
                ##########################################
                #mx_fig.suptitle('TROPoe '+version+' MIXR',x = plt_param.DLmeas_x_title,\
                mx_fig.suptitle(' MIXR',x = plt_param.DLmeas_x_title,\
                                y = plt_param.DLmeas_y_title2,horizontalalignment='right',\
                        fontsize=plt_param.TITLEfs);
                ###################################################################
                if (version == 'ASSIST'):
                    mx = mx_ax.pcolormesh(MEAStme_mesh,MEAShgt_mesh,MEASmx,vmin = MIXRmin,\
                                            vmax = MIXRmax,shading='nearest', cmap = cmap_mx,alpha=0.5) 
                    mx2 = mx_ax.pcolormesh(MEAStme_mesh,MEAShgt_mesh,MEASmx_nan,vmin = MIXRmin,\
                                            vmax = MIXRmax,shading='nearest', cmap = cmap_mx) 
                else:
                    mx2 = mx_ax.pcolormesh(MEAStme_mesh,MEAShgt_mesh,MEASmx,vmin = MIXRmin,\
                                            vmax = MIXRmax,shading='nearest', cmap = cmap_mx) 
                mx_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5);
                ###################################################################
                cbar = mx_fig.colorbar(mx2,ax=mx_ax,ticks = MIXRtick,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Mixing Ratio, g $\mathrm{kg}^{-1}$', fontsize=plt_param.LABELfs);
                ###################################################################
                # overlay cloud base height 
                mx_ax.plot(MEAStme_mesh[0],MEAScbh_nondef,'o',markeredgewidth=plt_param.PLTmewidth,\
                    markerfacecolor=plt_param.PLTmfcolor,markeredgecolor=plt_param.PLTmecolor,\
                    markersize=plt_param.PLTmsize,label=GMLcbhlabel)
                mx_ax.plot(MEAStme_mesh[0],MEAScbh_def,'o',markeredgewidth=plt_param.PLTmewidth,\
                    markerfacecolor=plt_param.PLTmfcolor_miss,markeredgecolor=plt_param.PLTmecolor,\
                    markersize=plt_param.PLTmsize,label='Cloud Base Height (default value, true value unknown)')
                mx_ax.legend(loc=8,fontsize=plt_param.LEGfs_bottom,bbox_to_anchor = plt_param.LEGpad_extra_bottom,\
                             handletextpad=.2, markerscale = 1, ncol = 2,frameon = False);
                # overlay contours
                tc = mx_ax.contour(MEAStme_mesh,MEAShgt_mesh,MEASmx,levels=MIXRclabels,colors=ccolor,\
                    linewidths=plt_param.PLTcwidth,linestyles='solid')
                mx_ax.clabel(tc,colors=ccolor,inline=True,inline_spacing=-1,fontsize=plt_param.PLTclsize)
                ###################################################################
                mx_ax.set_ylim(HGTmin,HGTmax); mx_ax.set_yticks(HGTtick);mx_ax.set_yticklabels(HGTtick_label);
                mx_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                # set figure y-label
                mx_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                mx_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                ###################################################################
                #insert data date as title
                T = mx_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                T.set_position(plt_param.Tloc);
                ###############################
    
                ###################################################################
                # Plot the Model Mixing Ratio Time Histories
                ##########################################
                mx_mdl = mdl_ax.imshow(MEASmx*np.nan,extent=DLgrid,vmin=MIXRmin,vmax=MIXRmax,\
                                     origin='lower',aspect='auto', cmap = cmap_mx);   
                mdl_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                cbar = mx_fig.colorbar(mx_mdl,ax=mdl_ax,ticks = MIXRtick,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Mixing Ratio, g $\mathrm{kg}^{-1}$', fontsize=plt_param.LABELfs);
                ###################################################################
                mdl_ax.set_ylim(HGTmin,HGTmax); mdl_ax.set_yticks(HGTtick);mdl_ax.set_yticklabels(HGTtick_label);
                mdl_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                ###################################################################
                # set figure y-label
                mdl_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                mdl_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                ###################################################################
                # Denote that Model Data was Not Available
                mdl_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([HGTmin,MAXhgt]),MDLname+' Data Not Available',\
                         horizontalalignment='center',verticalalignment='bottom',color='red',\
                         fontweight='semibold',fontsize=plt_param.TITLEfs);
                ###################################################################
                M = mdl_ax.set_title(MDLname + ' MIXR'+MDLstr,\
                            {'fontsize': plt_param.TITLEfs},pad=1,horizontalalignment='right');  
                M.set_position(plt_param.SPTLt_loc);
                
                """ ###############################################################
                Define and Plot the Model Mixing Ratio Biases
                ############################################################### """
                # Plot the Model Mixing Ratio Bias
                ###################################################################
                norm=colors.BoundaryNorm(boundaries=MIXRb_bounds,ncolors=256)
                #bias_mx = bias_ax.imshow(MEASmx*np.nan,extent=DLgrid,vmin=MIXRb_min,vmax=MIXRb_max,\
                #                     origin='lower',aspect='auto', cmap = cmap_bias_mx);   
                bias_mx = bias_ax.imshow(MEASmx*np.nan,extent=DLgrid,\
                                     origin='lower',aspect='auto', cmap = cmap_bias_mx, norm=norm);   
                bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                cbar = mx_fig.colorbar(bias_mx,ax=bias_ax,ticks = MIXRb_tick,pad=plt_param.Cbar_pad,\
                    aspect = plt_param.CBaspect,drawedges=True)
                cbar.ax.minorticks_off()
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Mix Ratio, g $\mathrm{kg}^{-1}$', fontsize=plt_param.LABELfs);
                ###################################################################
                bias_ax.set_ylim(HGTmin,HGTmax); bias_ax.set_yticks(HGTtick);bias_ax.set_yticklabels(HGTtick_label);
                bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                ###################################################################
                # set figure y-label
                bias_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                ###################################################################
                bias_ax.set_xticks(Xtick); bias_ax.set_xticklabels(Xtick_lab,fontsize=plt_param.TICKfs);
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=UserWarning)
                    bias_ax.set_xlim(Xmin,Xmax)
                ###################################################################
                bias_ax.set_xlabel('Hour, UTC',{'fontsize': plt_param.LABELfs});
                bias_ax.xaxis.set_label_coords(plt_param.Xlab_pad_tri[0],plt_param.Xlab_pad_tri[1]);
                ###################################################################
                # Denote that Model Data was Not Available
                bias_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([HGTmin,MAXhgt]),MDLname+' Data Not Available',\
                         horizontalalignment='center',verticalalignment='bottom',color='red',\
                         fontweight='semibold',fontsize=plt_param.TITLEfs);
                ###################################################################
                B = bias_ax.set_title('MIXR Error (Mdl - Obs)',{'fontsize': plt_param.TITLEfs},pad=1,\
                             horizontalalignment='right');  
                B.set_position(plt_param.SPTLt_loc);
                ###################################################################
                mx_fig.subplots_adjust(left = 0.1,right = 0.96,bottom = 0.125, top = 0.95)
                ###################################################################
                bias_ax.annotate(Date_start_label,xy = plt_param.DATEpad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                bias_ax.annotate(Date_end_label,xy = plt_param.DATE_end_pad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                #######################################################################
                # Add some extra info
                bias_ax.annotate(MDLgridinfo,xy = plt_param.GRIDpad_left_tri,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'left',\
                                fontsize = plt_param.TITLEfs)
                bias_ax.annotate(ATTRstr,xy = plt_param.ATTRpad_ctr_dual,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'center',\
                                fontsize = plt_param.TITLEfs)
                #######################################################################
                # add a logo
                imax = mx_fig.add_axes(plt_param.LOGOaxcoords)
                imax.set_axis_off()
                imax.imshow(logo, aspect="equal")
                ###################################################################          
                #Save the Developed Figure
                FILout = os.path.join(IMGpth_model,FILpfx+MDLname+'_MIXR_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                mx_fig.savefig(FILout,dpi=300); plt.close('all');
    
                """ ###############################################################
                Potential Temperature (instruct to share x-axis)
                ############################################################### """
                theta_fig, (theta_ax,mdl_ax,bias_ax) = plt.subplots(3,1,sharex = True,figsize=(plt_param.FIGwdth,plt_param.FIGwdth),\
                         height_ratios=[2,2,1],\
                         facecolor = 'w',edgecolor ='w',clear='True')
                ###################################################################
                # Plot the Lidar Potential Temperature Time Histories
                ##########################################
                #theta_fig.suptitle('TROPoe '+version+' THETA',x = plt_param.DLmeas_x_title,\
                theta_fig.suptitle('THETA',x = plt_param.DLmeas_x_title,\
                                y = plt_param.DLmeas_y_title2,horizontalalignment='right',\
                        fontsize=plt_param.TITLEfs);
                ###################################################################
                if (version == 'ASSIST'):
                    theta = theta_ax.pcolormesh(MEAStme_mesh,MEAShgt_mesh,MEAStheta,vmin = THETAmin,\
                                            vmax = THETAmax,shading='nearest', cmap = cmap_theta,alpha=0.5) 
                    theta2 = theta_ax.pcolormesh(MEAStme_mesh,MEAShgt_mesh,MEAStheta_nan,vmin = THETAmin,\
                                            vmax = THETAmax,shading='nearest', cmap = cmap_theta) 
                else:
                    theta2 = theta_ax.pcolormesh(MEAStme_mesh,MEAShgt_mesh,MEAStheta,vmin = THETAmin,\
                                            vmax = THETAmax,shading='nearest', cmap = cmap_theta) 
                theta_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5);
                ###################################################################
                cbar = theta_fig.colorbar(theta2,ax=theta_ax,ticks = THETAtick,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Potential Temperature, K', fontsize=plt_param.LABELfs);
                ###################################################################
                # overlay cloud base height 
                theta_ax.plot(MEAStme_mesh[0],MEAScbh_nondef,'o',markeredgewidth=plt_param.PLTmewidth,\
                    markerfacecolor=plt_param.PLTmfcolor,markeredgecolor=plt_param.PLTmecolor,\
                    markersize=plt_param.PLTmsize,label=GMLcbhlabel)
                theta_ax.plot(MEAStme_mesh[0],MEAScbh_def,'o',markeredgewidth=plt_param.PLTmewidth,\
                    markerfacecolor=plt_param.PLTmfcolor_miss,markeredgecolor=plt_param.PLTmecolor,\
                    markersize=plt_param.PLTmsize,label='Cloud Base Height (default value, true value unknown)')
                theta_ax.legend(loc=8,fontsize=plt_param.LEGfs_bottom,bbox_to_anchor = plt_param.LEGpad_extra_bottom,\
                             handletextpad=.2, markerscale = 1, ncol = 2,frameon = False);
                # overlay contours
                tc = theta_ax.contour(MEAStme_mesh,MEAShgt_mesh,MEAStheta,levels=THETAclabels,colors=ccolor,\
                    linewidths=plt_param.PLTcwidth,linestyles='solid')
                theta_ax.clabel(tc,colors=ccolor,inline=True,inline_spacing=-1,fontsize=plt_param.PLTclsize)
                ###################################################################
                theta_ax.set_ylim(HGTmin,HGTmax); theta_ax.set_yticks(HGTtick);theta_ax.set_yticklabels(HGTtick_label);
                theta_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                # set figure y-label
                theta_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                theta_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                ###################################################################
                #insert data date as title
                T = theta_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                T.set_position(plt_param.Tloc);
                ###############################
    
                ###################################################################
                # Plot the Model Potential Temperature Time Histories
                ##########################################
                theta_mdl = mdl_ax.imshow(MEAStheta*np.nan,extent=DLgrid,vmin=THETAmin,vmax=THETAmax,\
                                     origin='lower',aspect='auto', cmap = cmap_theta);   
                mdl_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                cbar = theta_fig.colorbar(theta_mdl,ax=mdl_ax,ticks = THETAtick,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Potential Temperature, K', fontsize=plt_param.LABELfs);
                ###################################################################
                mdl_ax.set_ylim(HGTmin,HGTmax); mdl_ax.set_yticks(HGTtick);mdl_ax.set_yticklabels(HGTtick_label);
                mdl_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                ###################################################################
                # set figure y-label
                mdl_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                mdl_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                ###################################################################
                # Denote that Model Data was Not Available
                mdl_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([HGTmin,MAXhgt]),MDLname+' Data Not Available',\
                         horizontalalignment='center',verticalalignment='bottom',color='red',\
                         fontweight='semibold',fontsize=plt_param.TITLEfs);
                ###################################################################
                M = mdl_ax.set_title(MDLname + ' THETA'+MDLstr,\
                            {'fontsize': plt_param.TITLEfs},pad=1,horizontalalignment='right');  
                M.set_position(plt_param.SPTLt_loc);
                
                """ ###############################################################
                Define and Plot the Model Potential Temperature Biases
                ############################################################### """
                # Plot the Model Potential Temperature Bias
                ###################################################################
                norm=colors.BoundaryNorm(boundaries=THETAb_bounds,ncolors=256)
                #bias_theta = bias_ax.imshow(MEAStheta*np.nan,extent=DLgrid,vmin=THETAb_min,vmax=THETAb_max,\
                #                     origin='lower',aspect='auto', cmap = cmap_bias_theta);   
                bias_theta = bias_ax.imshow(MEAStheta*np.nan,extent=DLgrid,\
                                     origin='lower',aspect='auto', cmap = cmap_bias_theta, norm=norm);   
                bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                cbar = theta_fig.colorbar(bias_theta,ax=bias_ax,ticks = THETAb_tick,pad=plt_param.Cbar_pad,\
                    aspect = plt_param.CBaspect,drawedges=True)
                cbar.ax.minorticks_off()
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Potential Temp, K', fontsize=plt_param.LABELfs);
                ###################################################################
                bias_ax.set_ylim(HGTmin,HGTmax); bias_ax.set_yticks(HGTtick);bias_ax.set_yticklabels(HGTtick_label);
                bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
                ###################################################################
                # set figure y-label
                bias_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                ###################################################################
                bias_ax.set_xticks(Xtick); bias_ax.set_xticklabels(Xtick_lab,fontsize=plt_param.TICKfs);
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=UserWarning)
                    bias_ax.set_xlim(Xmin,Xmax)
                ###################################################################
                bias_ax.set_xlabel('Hour, UTC',{'fontsize': plt_param.LABELfs});
                bias_ax.xaxis.set_label_coords(plt_param.Xlab_pad_tri[0],plt_param.Xlab_pad_tri[1]);
                ###################################################################
                # Denote that Model Data was Not Available
                bias_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([HGTmin,MAXhgt]),MDLname+' Data Not Available',\
                         horizontalalignment='center',verticalalignment='bottom',color='red',\
                         fontweight='semibold',fontsize=plt_param.TITLEfs);
                ###################################################################
                B = bias_ax.set_title('THETA Error (Mdl - Obs)',{'fontsize': plt_param.TITLEfs},pad=1,\
                             horizontalalignment='right');  
                B.set_position(plt_param.SPTLt_loc);
                ###################################################################
                theta_fig.subplots_adjust(left = 0.1,right = 0.96,bottom = 0.125, top = 0.95)
                ###################################################################
                bias_ax.annotate(Date_start_label,xy = plt_param.DATEpad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                bias_ax.annotate(Date_end_label,xy = plt_param.DATE_end_pad_sptl,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'right',\
                                fontsize = plt_param.TITLEfs)
                #######################################################################
                # Add some extra info
                bias_ax.annotate(MDLgridinfo,xy = plt_param.GRIDpad_left_tri,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'left',\
                                fontsize = plt_param.TITLEfs)
                bias_ax.annotate(ATTRstr,xy = plt_param.ATTRpad_ctr_dual,xycoords = 'figure fraction',\
                                verticalalignment = 'bottom',horizontalalignment = 'center',\
                                fontsize = plt_param.TITLEfs)
                #######################################################################
                # add a logo
                imax = theta_fig.add_axes(plt_param.LOGOaxcoords)
                imax.set_axis_off()
                imax.imshow(logo, aspect="equal")
                ###################################################################          
                #Save the Developed Figure
                FILout = os.path.join(IMGpth_model,FILpfx+MDLname+'_THETA_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                theta_fig.savefig(FILout,dpi=300); plt.close('all');
