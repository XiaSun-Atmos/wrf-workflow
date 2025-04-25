#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This Function Plots for the Case Wherein No Measurement Data Exists
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
######################################
from time_tuple_cat import time_cat
from atmo_temperature import calc_theta
from atmo_temperature import degC_to_K
###################################################
from atmo_ylim_plot import atmo_ylim
#####################################################
from cmaps import get_var_cmap
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
MAXhgt = 3000; 
ccolor = 'black'

def tropoe_model(MDLin,OUTinfo,version='ASSIST'):
    DATE = OUTinfo[3];
    basedate=datetime.strptime(DATE,'%Y%m%d')
    ###########################################################################
    IMGpth_meas = os.path.join(OUTinfo[0], 'no_meas',DATE);
    if not os.path.exists(IMGpth_meas):
        os.makedirs(IMGpth_meas)
    IMGpth_ext = os.path.join(OUTinfo[0],'ext', DATE);
    if not os.path.exists(IMGpth_ext):
        os.makedirs(IMGpth_ext)
    IMGpth_none = os.path.join(OUTinfo[0],'no_data', DATE);
    if not os.path.exists(IMGpth_none):
        os.makedirs(IMGpth_none)
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
    tropoe_Aplot(TROPoe_plt,HRRRw_plt,FILout)
    ###########################################################################
    MDLin ::  HRRRv4_plt = (HRRRv4_ini,HRRRv4_loc,HRRRv4_xy,HRRRv4_trp)
    OUTinfo :: FILout = (IMGout,MDLname,Sid,doi)
    ####################################################################### """
    ###########################################################################
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
                    ESRLtpmx_o = (To,To_sigma,Tsfc_o,Po,Psfc_o,MXrat_o,MXrat_sfc_o,Prate_o)
                    ESRLtpmx_int = (Tint,Tsfc_int,Pint,Psfc_int,MXrat_int,MXrat_sfc_int,Prate_int)
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

        ###########################################################################
        ###########################################################################
        """ Measured Wind Information Not Available -- Do Not Unpack """
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
                #MDLws = MDLws[:,Gref,:]; MDLwd = MDLwd[:,Gref,:];
                MDLt = MDLt[:,Gref,:]; MDLp = MDLp[:,Gref,:]; MDLmx = MDLmx[:,Gref,:];
                MDLtheta = calc_theta(degC_to_K(MDLt),MDLp)
                ###################################################################
                # Account for the way python interprets matrices (i.e. row x column)
                # ---> transpose inputted model data
                MDLt = np.transpose(MDLt); MDLmx = np.transpose(MDLmx);
                MDLtheta = np.transpose(MDLtheta)
                # now correctly (hgt x time) and consistent with the measured data
                
                ###################################################################
                """ Define Relevant Plotting Specifics 
                --> Specifics will vary depending on whether an extended forecast 
                --> is being produced
                ############################################################### """
                # define maximum and minimum measurement height 
                HGTmin,HGTmax,HGTtick = atmo_ylim((np.asarray([0,MAXhgt])),plt_param.SHGTrnd_bse,plt_param.TICKint);
                HGTtick_label = [str(ht/1000) for ht in HGTtick]
                ###################################################################
                # Define the Max Height Ind
                MDLhgt_max = np.where(MDLhgt > MAXhgt)[0]
                
                ###################################################################
                if len(MDLtime) > FXlen:
                    """ Define Bounds for Extended Forecast Period """
                    ###############################################################
                    # Adjust Model Variables #
                    MDLt_ext = MDLt; MDLmx_ext = MDLmx;
                    MDLtheta_ext = MDLtheta;
                    MDLtime_ext = MDLtime;
                    ###############################################################
                    """ Perform Model Averaging for Plotting and Statistical Analyses """
                    MDLtme_ext_mesh, MDLhgt_ext_mesh = np.meshgrid(MDLtime_ext,MDLhgt); #define reference meshgrid
                    MDLgrid_ext = (MDLtime_ext[0],MDLtime_ext[-1],MDLhgt[0],MDLhgt[-1])
                    ###############################################################
                    ext_mdl_oi = np.logical_and(MDLhgt_ext_mesh >= HGTmin, MDLhgt_ext_mesh <= HGTmax)
                    ###############################################################
                    # Set ranges based on the data instead of being hardcoded
                    #TEMPmin_ext,TEMPmax_ext,TEMPtick_ext = atmo_ylim((MDLt[ext_mdl_oi],\
                    #                MDLt_ext[ext_mdl_oi]),plt_param.TEMPrnd_bse,plt_param.TICKint)    
                    #TEMPclabels_ext=np.arange(TEMPmin_ext,TEMPmax_ext,plt_param.TEMPcint)
                    #MIXRmin_ext,MIXRmax_ext,MIXRtick_ext = atmo_ylim((MDLt[ext_mdl_oi],\
                    #                MDLt_ext[ext_mdl_oi]),plt_param.MIXRrnd_bse,plt_param.TICKint)    
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
    
                    """ Define Bounds for the Model-Specific Standard Forecast Period """
                    ###############################################################
                    # Adjust Model Variables #
                    MDLt = MDLt[:,0:FXlen]; MDLmx = MDLmx[:,0:FXlen];
                    MDLtheta = MDLtheta[:,0:FXlen]; 
                    MDLtime = MDLtime[0:FXlen];
                    ###############################################################
                    """ Perform Model Averaging for Plotting and Statistical Analyses """
                    MDLtme_mesh, MDLhgt_mesh = np.meshgrid(MDLtime,MDLhgt); #define reference meshgrid
                    MDLgrid = (MDLtime[0],MDLtime[-1],MDLhgt[0],MDLhgt[-1])
                    ###############################################################
                    mdl_oi = np.logical_and(MDLhgt_mesh >= HGTmin, MDLhgt_mesh <= HGTmax)
                    ###############################################################
                    # Set ranges based on the data instead of being hardcoded
                    #TEMPmin,TEMPmax,TEMPtick = atmo_ylim((MDLt[mdl_oi],MDLt[mdl_oi]),\
                    #                        plt_param.TEMPrnd_bse,plt_param.TICKint)    
                    #TEMPclabels=np.arange(TEMPmin,TEMPmax,plt_param.TEMPcint)
                    #MIXRmin,MIXRmax,MIXRtick = atmo_ylim((MDLt[mdl_oi],MDLt[mdl_oi]),\
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
                    MDLgrid = (MDLtime[0],MDLtime[-1],MDLhgt[0],MDLhgt[-1])
                    ###############################################################
                    mdl_oi = np.logical_and(MDLhgt_mesh >= HGTmin, MDLhgt_mesh <= HGTmax)
                    ###############################################################
                    # Set ranges based on the data instead of being hardcoded
                    #TEMPmin,TEMPmax,TEMPtick = atmo_ylim((MDLt[mdl_oi],MDLt[mdl_oi]),\
                    #                        plt_param.TEMPrnd_bse,plt_param.TICKint)    
                    #TEMPclabels=np.arange(TEMPmin,TEMPmax,plt_param.TEMPcint)
                    #MIXRmin,MIXRmax,MIXRtick = atmo_ylim((MDLt[mdl_oi],MDLt[mdl_oi]),\
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
                ###################################################################
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=UserWarning)
                    temp = temp_ax.imshow(MDLt*np.nan,extent=MDLgrid,vmin=TEMPmin,vmax=TEMPmax,\
                                       origin='lower',aspect='auto', cmap = cmap_temp);   
                temp_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5);
                ###################################################################
                cbar = temp_fig.colorbar(temp,ax=temp_ax,ticks = TEMPtick,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Temperature, $\mathrm{deg}^{\circ}$ C', fontsize=plt_param.LABELfs);
                ###################################################################
                temp_ax.set_ylim(HGTmin,HGTmax); temp_ax.set_yticks(HGTtick);temp_ax.set_yticklabels(HGTtick_label);
                temp_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                # set figure y-label
                temp_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                temp_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                ###################################################################
                # Denote that Measurement Data was Not Available
                temp_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([HGTmin,MAXhgt]),'Measurement Data Not Available',\
                         horizontalalignment='center',verticalalignment='bottom',color='red',\
                         fontweight='semibold',fontsize=plt_param.TITLEfs);
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
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=UserWarning)
                    norm=colors.BoundaryNorm(boundaries=TEMPb_bounds,ncolors=256)
                    #bias_temp = bias_ax.imshow(MDLt*np.nan,extent=MDLgrid,vmin=TEMPb_min,vmax=TEMPb_max,\
                    #                 origin='lower',aspect='auto', cmap = cmap_bias_temp);   
                    bias_temp = bias_ax.imshow(MDLt*np.nan,extent=MDLgrid,\
                                     origin='lower',aspect='auto', cmap = cmap_bias_temp, norm=norm);   
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
                B = bias_ax.set_title('TEMP Error (Mdl - Obs)',{'fontsize': plt_param.TITLEfs},pad=1,\
                             horizontalalignment='right');  
                B.set_position(plt_param.SPTLt_loc);
                ###################################################################
                # Denote that Measurement Data was Not Available
                bias_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([HGTmin,MAXhgt]),'Measurement Data Not Available',\
                         horizontalalignment='center',verticalalignment='bottom',color='red',\
                         fontweight='semibold',fontsize=plt_param.TITLEfs);
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
                # Add some attribution
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
                FILout = os.path.join(IMGpth_meas,FILpfx+MDLname+'_TEMP_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
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
                ###################################################################
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=UserWarning)
                    mx = mx_ax.imshow(MDLmx*np.nan,extent=MDLgrid,vmin=MIXRmin,vmax=MIXRmax,\
                                       origin='lower',aspect='auto', cmap = cmap_mx);   
                mx_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5);
                ###################################################################
                cbar = mx_fig.colorbar(mx,ax=mx_ax,ticks = MIXRtick,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Mixing Ratio, g $\mathrm{kg}^{-1}$', fontsize=plt_param.LABELfs);
                ###################################################################
                mx_ax.set_ylim(HGTmin,HGTmax); mx_ax.set_yticks(HGTtick);mx_ax.set_yticklabels(HGTtick_label);
                mx_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                # set figure y-label
                mx_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                mx_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                ###################################################################
                # Denote that Measurement Data was Not Available
                mx_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([HGTmin,MAXhgt]),'Measurement Data Not Available',\
                         horizontalalignment='center',verticalalignment='bottom',color='red',\
                         fontweight='semibold',fontsize=plt_param.TITLEfs);
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
                with warnings.catch_warnings():
                    norm=colors.BoundaryNorm(boundaries=MIXRb_bounds,ncolors=256)
                    warnings.simplefilter("ignore", category=UserWarning)
                    #bias_mx = bias_ax.imshow(MDLmx*np.nan,extent=MDLgrid,vmin=MIXRb_min,vmax=MIXRb_max,\
                    #                 origin='lower',aspect='auto', cmap = cmap_bias_mx);   
                    bias_mx = bias_ax.imshow(MDLmx*np.nan,extent=MDLgrid,\
                                     origin='lower',aspect='auto', cmap = cmap_bias_mx, norm=norm);   
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
                #bias_ax.set_xticks(Xtick); bias_ax.set_xticklabels(Xtick_lab,{'fontsize': plt_param.TICKfs});
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
                # Denote that Measurement Data was Not Available
                bias_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([HGTmin,MAXhgt]),'Measurement Data Not Available',\
                         horizontalalignment='center',verticalalignment='bottom',color='red',\
                         fontweight='semibold',fontsize=plt_param.TITLEfs);
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
                # Add some attribution
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
                FILout = os.path.join(IMGpth_meas,FILpfx+MDLname+'_MIXR_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
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
                ###################################################################
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=UserWarning)
                    theta = theta_ax.imshow(MDLtheta*np.nan,extent=MDLgrid,vmin=THETAmin,vmax=THETAmax,\
                                       origin='lower',aspect='auto', cmap = cmap_theta);   
                theta_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5);
                ###################################################################
                cbar = theta_fig.colorbar(theta,ax=theta_ax,ticks = THETAtick,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                cbar.set_label('Potential Temperature, K', fontsize=plt_param.LABELfs);
                ###################################################################
                theta_ax.set_ylim(HGTmin,HGTmax); theta_ax.set_yticks(HGTtick);theta_ax.set_yticklabels(HGTtick_label);
                theta_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                # set figure y-label
                theta_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                theta_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                ###################################################################
                # Denote that Measurement Data was Not Available
                theta_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([HGTmin,MAXhgt]),'Measurement Data Not Available',\
                         horizontalalignment='center',verticalalignment='bottom',color='red',\
                         fontweight='semibold',fontsize=plt_param.TITLEfs);
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
                with warnings.catch_warnings():
                    norm=colors.BoundaryNorm(boundaries=THETAb_bounds,ncolors=256)
                    warnings.simplefilter("ignore", category=UserWarning)
                    #bias_theta = bias_ax.imshow(MDLtheta*np.nan,extent=MDLgrid,vmin=THETAb_min,vmax=THETAb_max,\
                    #                 origin='lower',aspect='auto', cmap = cmap_bias_theta);   
                    bias_theta = bias_ax.imshow(MDLtheta*np.nan,extent=MDLgrid,\
                                     origin='lower',aspect='auto', cmap = cmap_bias_theta, norm=norm);   
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
                #bias_ax.set_xticks(Xtick); bias_ax.set_xticklabels(Xtick_lab,{'fontsize': plt_param.TICKfs});
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
                # Denote that Measurement Data was Not Available
                bias_ax.text(np.nanmean([Xmin,Xmax]),np.nanmean([HGTmin,MAXhgt]),'Measurement Data Not Available',\
                         horizontalalignment='center',verticalalignment='bottom',color='red',\
                         fontweight='semibold',fontsize=plt_param.TITLEfs);
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
                # Add some attribution
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
                FILout = os.path.join(IMGpth_meas,FILpfx+MDLname+'_THETA_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
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
                    temp = temp_ax.imshow(MDLt_ext*np.nan,extent=MDLgrid_ext,vmin=TEMPmin_ext,vmax=TEMPmax_ext,\
                                         origin='lower',aspect='auto', cmap = cmap_temp);   
                    temp_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5);
                    ###############################################################
                    cbar = temp_fig.colorbar(temp,ax=temp_ax,ticks = TEMPtick_ext,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                    cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                    cbar.set_label('Temperature, $\mathrm{deg}^{\circ}$ C', fontsize=plt_param.LABELfs);
                    ###############################################################
                    temp_ax.set_ylim(HGTmin,HGTmax); temp_ax.set_yticks(HGTtick);temp_ax.set_yticklabels(HGTtick_label);
                    temp_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                    ###################################################################
                    # set figure y-label
                    temp_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                    temp_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                    ###############################################################
                    # Denote that Measurement Data was Not Available
                    temp_ax.text(np.nanmean([Xmin_ext,Xmax_ext]),np.nanmean([HGTmin,MAXhgt]),'Measurement Data Not Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
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
                    norm=colors.BoundaryNorm(boundaries=TEMPb_bounds,ncolors=256)
                    #bias_temp = bias_ax.imshow(MDLt_ext*np.nan,extent=MDLgrid_ext,vmin=TEMPb_min,vmax=TEMPb_max,\
                    #                 origin='lower',aspect='auto', cmap = cmap_bias_temp);   
                    bias_temp = bias_ax.imshow(MDLt_ext*np.nan,extent=MDLgrid_ext,\
                                     origin='lower',aspect='auto', cmap = cmap_bias_temp, norm=norm);   
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
                    # Denote that Measurement Data was Not Available
                    bias_ax.text(np.nanmean([Xmin_ext,Xmax_ext]),np.nanmean([HGTmin,MAXhgt]),'Measurement Data Not Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
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
                    # Add some attribution
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
                    mx = mx_ax.imshow(MDLmx_ext*np.nan,extent=MDLgrid_ext,vmin=MIXRmin_ext,vmax=MIXRmax_ext,\
                                         origin='lower',aspect='auto', cmap = cmap_mx);   
                    mx_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5);
                    ###############################################################
                    cbar = mx_fig.colorbar(mx,ax=mx_ax,ticks = MIXRtick_ext,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                    cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                    cbar.set_label('Mixing Ratio, g $\mathrm{kg}^{-1}$', fontsize=plt_param.LABELfs);
                    ###############################################################
                    mx_ax.set_ylim(HGTmin,HGTmax); mx_ax.set_yticks(HGTtick);mx_ax.set_yticklabels(HGTtick_label);
                    mx_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                    ###################################################################
                    # set figure y-label
                    mx_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                    mx_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                    ###############################################################
                    # Denote that Measurement Data was Not Available
                    mx_ax.text(np.nanmean([Xmin_ext,Xmax_ext]),np.nanmean([HGTmin,MAXhgt]),'Measurement Data Not Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
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
                    norm=colors.BoundaryNorm(boundaries=MIXRb_bounds,ncolors=256)
                    #bias_mx = bias_ax.imshow(MDLmx_ext*np.nan,extent=MDLgrid_ext,vmin=MIXRb_min,vmax=MIXRb_max,\
                    #                 origin='lower',aspect='auto', cmap = cmap_bias_mx);   
                    bias_mx = bias_ax.imshow(MDLmx_ext*np.nan,extent=MDLgrid_ext,\
                                     origin='lower',aspect='auto', cmap = cmap_bias_mx, norm=norm);   
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
                    # Denote that Measurement Data was Not Available
                    bias_ax.text(np.nanmean([Xmin_ext,Xmax_ext]),np.nanmean([HGTmin,MAXhgt]),'Measurement Data Not Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
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
                    # Add some attribution
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
                    mx_fig, (mx_ax,mdl_ax,bias_ax) = plt.subplots(3,1,sharex = True,figsize=(plt_param.FIGwdth,plt_param.FIGwdth),\
                             height_ratios=[2,2,1],\
                             facecolor = 'w',edgecolor ='w',clear='True')
                    ###############################################################
                    # Plot the Potential Temperature Time Histories
                    ##########################################
                    #mx_fig.suptitle('TROPoe '+version+' MIXR',x = plt_param.DLmeas_x_title,\
                    mx_fig.suptitle(' MIXR',x = plt_param.DLmeas_x_title,\
                                    y = plt_param.DLmeas_y_title2,horizontalalignment='right',\
                            fontsize=plt_param.TITLEfs);
                    ###############################################################
                    mx = mx_ax.imshow(MDLmx_ext*np.nan,extent=MDLgrid_ext,vmin=MIXRmin_ext,vmax=MIXRmax_ext,\
                                         origin='lower',aspect='auto', cmap = cmap_mx);   
                    mx_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5);
                    ###############################################################
                    cbar = mx_fig.colorbar(mx,ax=mx_ax,ticks = MIXRtick_ext,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                    cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                    cbar.set_label('Potential Temperature, K', fontsize=plt_param.LABELfs);
                    ###############################################################
                    mx_ax.set_ylim(HGTmin,HGTmax); mx_ax.set_yticks(HGTtick);mx_ax.set_yticklabels(HGTtick_label);
                    mx_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                    ###################################################################
                    # set figure y-label
                    mx_ax.set_ylabel('Height, km AGL',{'fontsize': plt_param.LABELfs});
                    mx_ax.yaxis.set_label_coords(plt_param.Ylab_pad_dl[0],plt_param.Ylab_pad_dl[1]);
                    ###############################################################
                    # Denote that Measurement Data was Not Available
                    mx_ax.text(np.nanmean([Xmin_ext,Xmax_ext]),np.nanmean([HGTmin,MAXhgt]),'Measurement Data Not Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
                    ###############################################################
                    #insert data date as title
                    T = mx_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                    T.set_position(plt_param.Tloc);
                    ###############################
        
                    ###############################################################
                    # Plot the Model Potential Temperature Time Histories
                    ##########################################
                    mx_mdl = mdl_ax.pcolormesh(MDLtme_ext_mesh,MDLhgt_ext_mesh,MDLmx_ext,\
                                               vmin = MIXRmin_ext,vmax = MIXRmax_ext,shading='gouraud',cmap = cmap_mx)   
                    mdl_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                    ###############################################################
                    cbar = mx_fig.colorbar(mx_mdl,ax=mdl_ax,ticks = MIXRtick_ext,pad=plt_param.Cbar_pad,aspect = plt_param.CBaspect2)
                    cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                    cbar.set_label('Potential Temperature, K', fontsize=plt_param.LABELfs);
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
                    Define and Plot the Model Potential Temperature Biases
                    ########################################################### """
                    norm=colors.BoundaryNorm(boundaries=MIXRb_bounds,ncolors=256)
                    #bias_mx = bias_ax.imshow(MDLmx_ext*np.nan,extent=MDLgrid_ext,vmin=MIXRb_min,vmax=MIXRb_max,\
                    #                 origin='lower',aspect='auto', cmap = cmap_bias_mx);   
                    bias_mx = bias_ax.imshow(MDLmx_ext*np.nan,extent=MDLgrid_ext,\
                                     origin='lower',aspect='auto', cmap = cmap_bias_mx, norm=norm);   
                    bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                    ###############################################################
                    cbar = mx_fig.colorbar(bias_mx,ax=bias_ax,ticks = MIXRb_tick,pad=plt_param.Cbar_pad,\
                        aspect = plt_param.CBaspect,drawedges=True)
                    cbar.ax.minorticks_off()
                    cbar.ax.tick_params(labelsize=plt_param.Ctck_fs)
                    cbar.set_label('Potential Temp, K', fontsize=plt_param.LABELfs);
                    ###############################################################
                    bias_ax.set_ylim(HGTmin,HGTmax); bias_ax.set_yticks(HGTtick);bias_ax.set_yticklabels(HGTtick_label);
                    bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs);
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
                    # Denote that Measurement Data was Not Available
                    bias_ax.text(np.nanmean([Xmin_ext,Xmax_ext]),np.nanmean([HGTmin,MAXhgt]),'Measurement Data Not Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
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
                    # Add some attribution
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
    
                    del MDLtime_ext
      
