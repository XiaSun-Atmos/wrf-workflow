#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This Function Plots for the Case where both Measurement and Model Data are
Available
    
###############################################################################
Created on Fri Oct 25 13:31:49 2019
@author: jduncan
###############################################################################
The Purpose of this Script is to Develop and Output Wind-Related Plots Based on
Measurements from the SFC Standard Meteorological Information Installed at 
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
import os,time
from datetime import datetime,timedelta
#####################################
from coord_distance import dist2coord
from atmo_meas_to_hour_end  import atmo_meas2hr_end
from wd_meas_to_hour_end import wd_meas2hr_end
from time_tuple_cat import time_cat
from atmo_tuple_cat import atmo_cat
from atmo_time_correction import atmo_time_correct
#####################################
from atmo_bias import atmo_bias
from wd_bias_flip import wd_bias
#####################################
from atmo_ylim_plot import atmo_ylim
from atmo_ylim_bias import atmo_blim
from wd_ylim_plot import wd_ylim
from wd_ylim_bias import wd_blim
from wd_disc import wd_plt
from avgmins_to_label import avgmins_to_label
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
AVGper = getattr(plt_param,'PSLmet_avgmins',60)
AVGlabel = avgmins_to_label(AVGper)
###############################################################################
# Define SFC Measurement Height 
SFCws_hgt = 10;
# per Jenni Kyrouac -- For the met systems: winds = 10m, pressure = 1m, temp/rh = 2m.

def psl_Wplot(MEASin,MDLin,OUTinfo):
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
    psl_Wplot(MEASin,MDLin,FILout)
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
            MDLloc = (SGPsite,SGPlon,SGPlat,Mlon,Mlat,Melev,Mdist)
            MDLxy = (Mfcst_tme,Mmeas_hgt)
            MDLwind = depends on interpolation choice (see below)
                    ESRLw_o = (WSo,WDo,WSsfc_o,WDsfc_o,Uo,Usfc_o,Vo,Vsfc_o,Wo)
                    ESRLw_int = (WSint,WDint,WDsfc_o,WDsfc_int,Uint,Usfc_int,Vint,Vsfc_int,Wint)   
                    ###########################################################################
                    ESRLw = (ESRLw_o,ESRLw_int)
    ########################################################################"""
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
    
    ###########################################################################
    """ Unpack Measured Wind Information 
    ####################################
    MEASw_plt = (MEASloc,MEAStme,MEASw_in)         
    ---> MEASloc = (MEASlon,MEASlat,MEASalt)
    ---> MEAStme = (MEASxy_minus,MEASxy,MEASxy_plus)
    ---> MEASw_in = (MEASwind_minus,MEASwind,MEASwind_plus)
    ###########################################################################
    MEASw = (MEASws,MEASwd_vec)
    ########################################################################"""
    # Define Relevant Indices for Reference in Concatenation #
    WSind_cat = 0; WSqc_ind_cat = 1;
    WDind_cat = 1; WDqc_ind_cat = 1;
    ###########################################################################
    
    #-> Combine Relevant Measurement Tuples 
    MEAStime = time_cat(MEASin[1][0],MEASin[1][1],MEASin[1][2],MEASin[1][3])
    ###########################################################################
    MEASw_in = MEASin[2]; 
    #####################
    MEASws = atmo_cat(MEASw_in[0],MEASw_in[1],MEASw_in[2],MEASw_in[3],WSind_cat)
    MEASws_qc = None
    MEASwd = atmo_cat(MEASw_in[0],MEASw_in[1],MEASw_in[2],MEASw_in[3],WDind_cat)
    MEASwd_qc = None

    ###################################################################
    # Perform Correction to Create Equally Spaced Data Arrays #
    ###################################################################
    MEAStime_new, MEASws = atmo_time_correct(MEAStime,MEASws);
    _, MEASwd = atmo_time_correct(MEAStime,MEASwd);
    MEAStime = MEAStime_new;

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
            # Determine Whether Model Data Exists for Said Initialization Time
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
                """ Define Reference Height and Time Indice for Model Data """
                MDLhgt = MDLxy[INIref][1];
                ################################
                HGTdel = abs(MDLhgt - SFCws_hgt)
                MDLhgt_ref = np.where((HGTdel - np.nanmin(HGTdel)) < 1e-5)[0][0]  # index of level
    
                """ Extract Relevant Model Information """
                MDLw_ini = MDLwind[INIref][INTPoi]
                ###################################################################
                # SFC -- Interested in Surface Information |(forecast_hour,sites,hgt)|
                MDLws = MDLw_ini[0];
                MDLwd = MDLw_ini[1];
                MDLws_sfc = MDLw_ini[2];
                MDLwd_sfc = MDLw_ini[3];
                ###################################################################
                #MDLws_sfc = MDLw_ini[2] #surface wind speed information   
                #MDLwd_sfc = MDLw_ini[3] #surface wind direction information 
                ###################################################################
                MDLws = MDLws[:,Gref,MDLhgt_ref];
                MDLwd = MDLwd[:,Gref,MDLhgt_ref];
                MDLws_sfc = MDLws_sfc[:,Gref];
                MDLwd_sfc = MDLwd_sfc[:,Gref];
                #MDLws = MDLws_sfc
                #MDLwd = MDLwd_sfc
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
                    MDLws_ext = MDLws; MDLwd_ext = MDLwd;
                    MDLtime_ext = MDLtime;
                    ###############################################################
                    # Perform Model Averaging for Plotting and Statistical Analsyes
                    MEASws_ext_mu = atmo_meas2hr_end(MEASws,MEASws_qc,MEAStime,MDLtime_ext, AVGmins=AVGper);
                    MEASwd_ext_mu = wd_meas2hr_end(MEASwd,MEASwd_qc,MEAStime,MDLtime_ext, AVGmins=AVGper);
                    ###############################################################
                    #WSmin_ext,WSmax_ext,WStick_ext = atmo_ylim((MEASws[np.logical_and(MEASws_qc == 0,ext_meas_oi)],\
                    WSmin_ext,WSmax_ext,WStick_ext = atmo_ylim((MEASws[ext_meas_oi],\
                                                           MDLws_ext),plt_param.WSrnd_bse,plt_param.TICKint);
                    #WDmin_ext,WDmax_ext,WDtick_ext = wd_ylim((MEASwd[np.logical_and(MEASwd_qc == 0,ext_meas_oi)],MDLwd_ext));
                    WDmin_ext,WDmax_ext,WDtick_ext = wd_ylim((MEASwd[ext_meas_oi],MDLwd_ext));
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
                    MDLws = MDLws[0:FXlen]; MDLwd = MDLwd[0:FXlen];
                    MDLtime = MDLtime[0:FXlen];
                    ###############################################################
                    # Perform Model Averaging for Plotting and Statistical Analsyes
                    MEASws_mu = atmo_meas2hr_end(MEASws,MEASws_qc,MEAStime,MDLtime, AVGmins=AVGper);
                    MEASwd_mu = wd_meas2hr_end(MEASwd,MEASwd_qc,MEAStime,MDLtime, AVGmins=AVGper);
                    ###############################################################
                    #WSmin,WSmax,WStick = atmo_ylim((MEASws[np.logical_and(MEASws_qc == 0,meas_oi)],\
                    WSmin,WSmax,WStick = atmo_ylim((MEASws[meas_oi],\
                                                           MDLws),plt_param.WSrnd_bse,plt_param.TICKint);
                    #WDmin,WDmax,WDtick = wd_ylim((MEASwd[np.logical_and(MEASwd_qc == 0,meas_oi)],MDLwd));
                    WDmin,WDmax,WDtick = wd_ylim((MEASwd[meas_oi],MDLwd));
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
                    MEASws_mu = atmo_meas2hr_end(MEASws,MEASws_qc,MEAStime,MDLtime, AVGmins=AVGper);
                    MEASwd_mu = wd_meas2hr_end(MEASwd,MEASwd_qc,MEAStime,MDLtime, AVGmins=AVGper);
                    # desire mean values at each desired forecast time (push forward)
                    ###############################################################
                    #WSmin,WSmax,WStick = atmo_ylim((MEASws[np.logical_and(MEASws_qc == 0,meas_oi)],\
                    WSmin,WSmax,WStick = atmo_ylim((MEASws[meas_oi],\
                                                           MDLws),plt_param.WSrnd_bse,plt_param.TICKint);
                    #WDmin,WDmax,WDtick = wd_ylim((MEASwd[np.logical_and(MEASwd_qc == 0,meas_oi)],MDLwd));
                    WDmin,WDmax,WDtick = wd_ylim((MEASwd[meas_oi],MDLwd));
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
                #ws_ax.plot(MEAStime[np.where(MEASws_qc == 0)],MEASws[np.where(MEASws_qc == 0)],color='black',\
                ws_ax.plot(MDLtime,MEASws_mu,'s',color='gray',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Obs (1Hr Avg)')
                ws_ax.plot(MEAStime,MEASws,color='black', linewidth=plt_param.PLTlwidth,label='Obs')
                ws_ax.plot(MDLtime,MDLws,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname+MDLstr)
                ws_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                handles, labels = ws_ax.get_legend_handles_labels()
                order = [1,0,2]
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
                #insert location as title
                T = ws_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                T.set_position(plt_param.Tloc);
                
                """ ###########################################################
                Define and Plot the Model Wind Speed Biases
                ############################################################"""
                # initate 0/0 line
                bias_ax.plot([Xmin,Xmax],[0,0],color='black',linewidth=1.5)
                ###################################################################
                # Define the Bias (Observations - Measurements)
                if (np.isnan(MEASws_mu).all()):
                    bias_ax.plot(MDLtime,MDLtime*np.nan,'^-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Mdl - Obs')
                else:
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
                if (np.isnan(MEASws_mu).all()):
                    WSb_min=-5;WSb_max=5;WSb_tick=np.arange(WSb_min,WSb_max+1,5)
                else:
                    WSb_min,WSb_max,WSb_tick = atmo_blim((WSbias),plt_param.WSbias_bse,plt_param.TICKint)
                bias_ax.set_ylim(WSb_min,WSb_max); bias_ax.set_yticks(WSb_tick);
                bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                bias_ax.set_ylabel('Model Error, m $\mathrm{s}^{-1}$',{'fontsize': plt_param.LABELfs});    
                bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                ###################################################################
                #input basic statistics
                if (np.isnan(MEASws_mu).all()):
                    bias_ax.text(np.nanmean([Xmin,Xmax]),0,'Measurement Data Not Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
                else:
                    RMSEws,BIASws = atmo_bias(WSbias)
                    #insert data statistics as title
                    BT = bias_ax.set_title('RMSE: ' + str(RMSEws)+' m '+r'$\mathrm{s}^{-1}$ | Mean Bias: '+\
                                      str(BIASws)+' m '+r'$\mathrm{s}^{-1}$',{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5) 
                    BT.set_position(plt_param.Bloc_left);
                ###################################################################
                #-> annotate measurement heights
                ws_ax.annotate('10 m AGL',xy=(0.13,0.895), xycoords='figure fraction',color='black',\
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
                #######################################################################
                # add a logo
                imax = ws_fig.add_axes(plt_param.metLOGOaxcoords)
                imax.set_axis_off()
                imax.imshow(logo, aspect="equal")
                ###################################################################
                #Save the Developed Figure
                FILout = os.path.join(IMGpth,'pslmet_'+MDLname+'_WS_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                ws_fig.savefig(FILout,dpi=300); plt.close('all');
                
                """ ###############################################################
                Wind Direction Plots (instruct to share x-axis)
                ############################################################### """
                wd_fig, (wd_ax,bias_ax) = plt.subplots(2,1,sharex = True, figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                         facecolor = 'w',edgecolor ='w',clear='True')
                ###################################################################
                # Plot the Wind Direction Time Histories
                ###################################################################
                wd_ax = wd_plt(MDLtime,MEASwd_mu,wd_ax,'s','gray','Obs (1Hr Avg)',marker_size=plt_param.PLTmsize)
                wd_ax = wd_plt(MEAStime,MEASwd,wd_ax,'.','black','Obs')
                wd_ax = wd_plt(MDLtime,MDLwd,wd_ax,'o-','red',MDLname+MDLstr,marker_size=plt_param.PLTmsize,line_width=plt_param.PLTlwidth)
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
                #insert location as title
                T = wd_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                T.set_position(plt_param.Tloc);
                
                """ ###############################################################
                # Plot the Model Wind Direction Biases
                ############################################################### """
                # initate 0/0 line
                bias_ax.plot([Xmin,Xmax],[0,0],color='black',linewidth=1.5)
                ###################################################################
                # Define the Bias (Observations - Measurements)
                if (np.isnan(MEASwd_mu).all()):
                    bias_ax.plot(MDLtime,MDLtime*np.nan,'^-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Mdl - Obs')
                else:
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
                if (np.isnan(MEASwd_mu).all()):
                    WDb_min=-180;WDb_max=180;WDb_tick = np.arange(WDb_min,WDb_max+1,60);
                else:
                    WDb_min,WDb_max,WDb_tick = wd_blim((WDbias));
                bias_ax.set_ylim(WDb_min,WDb_max); bias_ax.set_yticks(WDb_tick);
                bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                bias_ax.set_ylabel('Model Error, deg ${}^{\circ}$',{'fontsize': plt_param.LABELfs});    
                bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                ###################################################################
                #insert data statistics as title
                if (np.isnan(MEASwd_mu).all()):
                    bias_ax.text(np.nanmean([Xmin,Xmax]),0,'Measurement Data Not Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
                else:
                    RMSEws,BIASws = atmo_bias(WSbias)
                    #insert data statistics as title
                    BT = bias_ax.set_title('RMSE: ' + str(RMSEwd)+' deg '+r'$\mathrm{}^{\circ}$ | Mean Bias: '+\
                                      str(BIASwd)+' deg '+r'$\mathrm{}^{\circ}$',{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                    BT.set_position(plt_param.Bloc_left);
                ###################################################################
                #-> annotate measurement heights
                wd_ax.annotate('10 m AGL',xy=(0.13,0.895), xycoords='figure fraction',color='black',\
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
                #######################################################################
                # add a logo
                imax = wd_fig.add_axes(plt_param.metLOGOaxcoords)
                imax.set_axis_off()
                imax.imshow(logo, aspect="equal")
                ###################################################################
                #Save the Developed Figure
                FILout = os.path.join(IMGpth,'pslmet_'+MDLname+'_WD_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                wd_fig.savefig(FILout,dpi=300); plt.close('all')
                ###################################################################
                
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
                    #ws_ax.plot(MEAStime[np.where(MEASws_qc == 0)],MEASws[np.where(MEASws_qc == 0)],color='black',\
                    ws_ax.plot(MDLtime_ext,MEASws_ext_mu,'s',color='gray',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Obs (1Hr Avg)')
                    ws_ax.plot(MEAStime,MEASws,color='black', linewidth=plt_param.PLTlwidth,label='Obs')
                    ws_ax.plot(MDLtime_ext,MDLws_ext,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname+MDLstr)
                    ws_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                    ###############################################################
                    handles, labels = ws_ax.get_legend_handles_labels()
                    order = [1,0,2]
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
                    #insert location as title
                    T = ws_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                    T.set_position(plt_param.Tloc);
                    
                    """ ###########################################################
                    Define and Plot the Model Wind Speed Biases
                    ########################################################### """
                    # initate 0/0 line
                    bias_ax.plot([Xmin_ext,Xmax_ext],[0,0],color='black',linewidth=1.5)
                    ###############################################################
                    # Define the Bias (Observations - Measurements)
                    if (np.isnan(MEASws_ext_mu).all()):
                        bias_ax.plot(MDLtime_ext,MDLtime_ext*np.nan,'^-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Mdl - Obs')
                    else:
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
                    if (np.isnan(MEASws_ext_mu).all()):
                        WSb_min = -5;WSb_max = 5;WSb_tick = np.arange(WSb_min,WSb_max+1,5)
                    else:
                        WSb_min,WSb_max,WSb_tick = atmo_blim((WSbias),plt_param.WSbias_bse,plt_param.TICKint)
                    bias_ax.set_ylim(WSb_min,WSb_max); bias_ax.set_yticks(WSb_tick);
                    bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                    ###############################################################
                    bias_ax.set_ylabel('Model Error, m $\mathrm{s}^{-1}$',{'fontsize': plt_param.LABELfs});    
                    bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                    ###############################################################
                    #input basic statistics
                    if (np.isnan(MEASws_ext_mu).all()):
                        bias_ax.text(np.nanmean([Xmin_ext,Xmax_ext]),0,'Measurement Data Not Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
                    else:
                        RMSEws,BIASws = atmo_bias(WSbias)
                        #insert data statistics as title
                        BT = bias_ax.set_title('RMSE: ' + str(RMSEws)+' m '+r'$\mathrm{s}^{-1}$ | Mean Bias: '+\
                                          str(BIASws)+' m '+r'$\mathrm{s}^{-1}$',{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5) 
                        BT.set_position(plt_param.Bloc_left);
                    ###############################################################
                    #-> annotate measurement heights
                    ws_ax.annotate('10 m AGL',xy=(0.13,0.895), xycoords='figure fraction',color='black',\
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
                    #######################################################################
                    # add a logo
                    imax = ws_fig.add_axes(plt_param.metLOGOaxcoords)
                    imax.set_axis_off()
                    imax.imshow(logo, aspect="equal")
                    ###################################################################
                    #Save the Developed Figure
                    FILout = os.path.join(IMGpth_ext,'pslmet_'+MDLname+'_WS_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                    ws_fig.savefig(FILout,dpi=300); plt.close('all');
                    
                    """ ###########################################################
                    Wind Direction Plots (instruct to share x-axis)
                    ########################################################### """
                    wd_fig, (wd_ax,bias_ax) = plt.subplots(2,1,sharex = True, figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                             facecolor = 'w',edgecolor ='w',clear='True')
                    ###############################################################
                    # Plot the Wind Direction Time Histories
                    ###############################################################
                    wd_ax = wd_plt(MDLtime_ext,MEASwd_ext_mu,wd_ax,'s','gray','Obs (1Hr Avg)',marker_size=plt_param.PLTmsize)
                    wd_ax = wd_plt(MEAStime,MEASwd,wd_ax,'.','black','Obs')
                    wd_ax = wd_plt(MDLtime_ext,MDLwd_ext,wd_ax,'o-','red',MDLname+MDLstr,marker_size=plt_param.PLTmsize,line_width=plt_param.PLTlwidth)
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
                    #insert location as title
                    T = wd_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                    T.set_position(plt_param.Tloc);
                    
                    """ ###########################################################
                    # Plot the Model Wind Direction Biases
                    ########################################################### """
                    # initate 0/0 line
                    bias_ax.plot([Xmin_ext,Xmax_ext],[0,0],color='black',linewidth=1.5)
                    ###############################################################
                    # Define the Bias (Observations - Measurements)
                    if (np.isnan(MEASwd_ext_mu).all()):
                        bias_ax.plot(MDLtime_ext,MDLtime_ext*np.nan,'^-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Mdl - Obs')
                    else:
                        WDbias,RMSEwd,BIASwd = wd_bias(MEASwd_ext_mu,MDLwd_ext)
                        bias_ax.plot(MDLtime_ext,WDbias,'^-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Mdl - Obs')
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
                    #Define Bounds of Wind Speed - Model Bias to Determine Y-Bounds
                    if (np.isnan(MEASwd_ext_mu).all()):
                        WDb_min=-180;WDb_max=180;WDb_tick=np.arange(WDb_min, WDb_max+1,60);
                    else:
                        WDb_min,WDb_max,WDb_tick = wd_blim((WDbias));
                    bias_ax.set_ylim(WDb_min,WDb_max); bias_ax.set_yticks(WDb_tick);
                    bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                    ###############################################################
                    bias_ax.set_ylabel('Model Error, deg ${}^{\circ}$',{'fontsize': plt_param.LABELfs});    
                    bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                    ###############################################################
                    #insert data statistics as title
                    if (np.isnan(MEASwd_ext_mu).all()):
                        bias_ax.text(np.nanmean([Xmin_ext,Xmax_ext]),0,'Measurement Data Not Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
                    else:
                        BT = bias_ax.set_title('RMSE: ' + str(RMSEwd)+' deg '+r'$\mathrm{}^{\circ}$ | Mean Bias: '+\
                                      str(BIASwd)+' deg '+r'$\mathrm{}^{\circ}$', \
                                      {'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                        BT.set_position(plt_param.Bloc_left);
                    ###############################################################
                    #-> annotate measurement heights
                    wd_ax.annotate('10 m AGL',xy=(0.13,0.895), xycoords='figure fraction',color='black',\
                                     fontsize=plt_param.LEGfs,weight='bold')
                    ###############################################################
                    wd_fig.subplots_adjust(left = 0.125,right = 0.935,bottom = 0.155, top = 0.925)
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
                    #######################################################################
                    # add a logo
                    imax = wd_fig.add_axes(plt_param.metLOGOaxcoords)
                    imax.set_axis_off()
                    imax.imshow(logo, aspect="equal")
                    ###################################################################
                    #Save the Developed Figure
                    FILout = os.path.join(IMGpth_ext,'pslmet_'+MDLname+'_WD_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                    wd_fig.savefig(FILout,dpi=300); plt.close('all')
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
                MEASws_mu = atmo_meas2hr_end(MEASws,MEASws_qc,MEAStime,MDLtime, AVGmins=AVGper);
                MEASwd_mu = wd_meas2hr_end(MEASwd,MEASwd_qc,MEAStime,MDLtime, AVGmins=AVGper);
                
                """ Define Relevant Plotting Specifics """
                # Define Proper Indices to Determine Axes Parameters
                meas_oi = np.logical_and(MEAStime >= MDLtime[0], MEAStime <= MDLtime[-1]) 
                ###################################################################
                #WSmin,WSmax,WStick = atmo_ylim((MEASws[np.logical_and(MEASws_qc == 0,meas_oi)]),\
                WSmin,WSmax,WStick = atmo_ylim((MEASws[meas_oi]),\
                                               plt_param.WSrnd_bse,plt_param.TICKint)
                #WDmin,WDmax,WDtick = wd_ylim((MEASwd[np.logical_and(MEASwd_qc == 0,meas_oi)]));
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
                #ws_ax.plot(MEAStime[np.where(MEASws_qc == 0)],MEASws[np.where(MEASws_qc == 0)],color='black',\
                ws_ax.plot(MDLtime,MEASws_mu,'s',color='gray',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label='Obs (1Hr Avg)')
                ws_ax.plot(MEAStime,MEASws,color='black', linewidth=plt_param.PLTlwidth,label='Obs')
                ws_ax.plot(MDLtime,MDLtime*np.nan,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname+MDLstr)
                ws_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                handles, labels = ws_ax.get_legend_handles_labels()
                order = [1,0,2]
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
                #insert location as title
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
                bias_ax.set_ylabel('Model Error, m $\mathrm{s}^{-1}$',{'fontsize': plt_param.LABELfs});    
                bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                ###################################################################
                #insert data statistics as title
                BT = bias_ax.set_title('RMSE: N/A',{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                BT.set_position(plt_param.Bloc_left);
                ###################################################################
                #-> annotate measurement heights
                ws_ax.annotate('10 m AGL',xy=(0.13,0.895), xycoords='figure fraction',color='black',\
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
                #######################################################################
                # add a logo
                imax = ws_fig.add_axes(plt_param.metLOGOaxcoords)
                imax.set_axis_off()
                imax.imshow(logo, aspect="equal")
                ###################################################################
                #Save the Developed Figure
                FILout = os.path.join(IMGpth_model,'pslmet_'+MDLname+'_WS_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                ws_fig.savefig(FILout,dpi=300); plt.close('all');
                        
                """ ###############################################################
                Wind Direction Plots (instruct to share x-axis)
                ############################################################### """
                wd_fig, (wd_ax,bias_ax) = plt.subplots(2,1,sharex = True, figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),\
                         facecolor = 'w',edgecolor ='w',clear='True')
                ###################################################################
                # Plot the Wind Direction Time Histories
                ###################################################################
                #MEASwd[MEASwd_qc != 0] = np.nan;
                wd_ax = wd_plt(MDLtime,MEASwd_mu,wd_ax,'s','gray','Obs (1Hr Avg)',marker_size=plt_param.PLTmsize)
                wd_ax = wd_plt(MEAStime,MEASwd,wd_ax,'.','black','Obs')
                wd_ax.plot(MDLtime,MDLtime*np.nan,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname+MDLstr)
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
                #insert location as title
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
                bias_ax.set_ylabel('Model Error, deg ${}^{\circ}$',{'fontsize': plt_param.LABELfs});    
                bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                ###################################################################
                #insert data statistics as title
                BT = bias_ax.set_title('RMSE: N/A',{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                BT.set_position(plt_param.Bloc_left);
                ###################################################################
                #-> annotate measurement heights
                wd_ax.annotate('10 m AGL',xy=(0.13,0.895), xycoords='figure fraction',color='black',\
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
                #######################################################################
                # add a logo
                imax = wd_fig.add_axes(plt_param.metLOGOaxcoords)
                imax.set_axis_off()
                imax.imshow(logo, aspect="equal")
                ###################################################################
                #Save the Developed Figure
                FILout = os.path.join(IMGpth_model,'pslmet_'+MDLname+'_WD_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                wd_fig.savefig(FILout,dpi=300); plt.close('all')
                ###################################################################
