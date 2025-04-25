#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This Function Plots for the Case where both Measurement and Model Data are 
Available

###############################################################################
Created on Fri Oct 25 13:31:49 2019
@author: jduncan
###############################################################################
The Purpose of this Script is to Develop and Output Temp/Pressure/Mixing Ratio
Plots Based on Measurements from the SFC Standard Meteorological Information 
Installed at Various Sites
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

This Function Will No Longer Output Plots of Relative Humidity -- Instead Plots
of the Mixing Ratio will be Produced
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as image
import os,time,warnings
from datetime import datetime,timedelta
#####################################
from coord_distance import dist2coord
from atmo_meas_to_hour_end import atmo_meas2hr_end
from atmo_meas_to_hour_upto import atmo_meas2hr_upto
from time_tuple_cat import time_cat
from atmo_tuple_cat import atmo_cat
from atmo_tuple_cat import atmo_cat
from atmo_time_correction import atmo_time_correct
from atmo_moisture import calc_dewpoint_from_RH
from atmo_moisture import calc_dewpoint_from_MR
from atmo_moisture import calc_mixingratio
from atmo_moisture import calc_rh_from_MR
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
AVGper = getattr(plt_param,'PSLmet_avgmins',60)
AVGlabel = avgmins_to_label(AVGper)
###############################################################################

def psl_RADmodel(MDLin,OUTinfo):
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
    FILpfx='pslmet_'
    
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
    psl_RADmodel(HRRRtrp_plt,FILout)
    #############################################
    OUTinfo :: FILout = (IMGout,MDLname,Sid,doi)
    #############################################"""

    ###########################################################################
    """ Unpack Modeled Temperature, Pressure, and Mixing Ratio Data
    ###############################################################
    MDLin :: 
        HRRRrad_plt = (ESRLhrrr_ini,ESRLhrrr_loc,ESRLhrrr_xy,ESRLhrrr_rad)
        RAPrad_plt = (ESRLrap_ini,ESRLrap_loc,ESRLrap_xy,ESRLrap_rad)
        #-> each embedded tuple has length of initialization files
            ###################################################################
            MDLloc = (MDLsite,MDLlon,MDLlat,Mlon,Mlat,Melev,Mdist)
            MDLxy = (Mfcst_tme,Mmeas_hgt)
            MDLrad = (ESRLrad_flux, ESRLrad_dd, ESRLrad_cs)
                ESRLrad_flux = (ESRLrad_o,ESRLrad_int)
                ###############################################################
                ESRLrad_o = (SWdsfc_o, LWdsfc_o,SWusfc_o, LWusfc_o)
                ESRLrad_int = (SWdsfc_int, LWdsfc_int,SWusfc_int, LWusfc_int)
    ########################################################################"""
    # Model Initialization Times #        
    MDLini = MDLin[0]; #outputs model initialization time for each ini hour
    MDLxy = MDLin[2]; #outputs forecast time and measurement height
    MDLsrad = MDLin[3] #outputs ini hour srad data (expand in loop)
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
    """ No measured RAD available
    ####################################"""

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
                MDLsrad_ini = MDLsrad[INIref]
                ##########################################
                # this INTPoi reference convention is correct
                MDLsrad_flux = MDLsrad_ini[0][INTPoi*2]
                MDLsrad_dd = MDLsrad_ini[1][INTPoi]
                ###################################################################
                # Extract Relevant Model Surface Radiation Information
                MDLsw_dwn = MDLsrad_flux[0]; MDLsw_upp = MDLsrad_flux[2];
                MDLlw_dwn = MDLsrad_flux[1]; MDLlw_upp = MDLsrad_flux[3];
    
                #direct/diffuse information
                MDLsw_dir = MDLsrad_dd[0]
                MDLsw_dif = MDLsrad_dd[1];
                ###################################################################
                # modify direct/diffuse radiation information to accomodate RAP
                # reported -999.0 values
                ###################################################################
                MDLsw_dir[MDLsw_dir == -999.0] = np.nan;
                MDLsw_dif[MDLsw_dif == -999.0] = np.nan;
                ###################################################################
                # Maintain Data Only From Relevant Grid Point
                MDLsw_dwn = MDLsw_dwn[:,Gref]; MDLsw_upp = MDLsw_upp[:,Gref];
                MDLlw_dwn = MDLlw_dwn[:,Gref]; MDLlw_upp = MDLlw_upp[:,Gref];

                MDLsw_dir = MDLsw_dir[:,Gref];
                MDLsw_dif = MDLsw_dif[:,Gref];

                MDLnet_rad = (MDLsw_dwn+MDLlw_dwn)-(MDLsw_upp+MDLlw_upp)
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
                    ###############################################################
                    # Adjust Model Variables #
                    MDLsw_dwn_ext = MDLsw_dwn; MDLnet_rad_ext = MDLnet_rad; 
                    MDLtime_ext = MDLtime;
                    ###############################################################
                    SRADmin_ext,SRADmax_ext,SRADtick_ext = atmo_ylim((MDLsw_dwn_ext),plt_param.SRADrnd_bse,plt_param.TICKint)
                    NRADmin_ext,NRADmax_ext,NRADtick_ext = atmo_ylim((MDLnet_rad_ext),plt_param.SRADrnd_bse,plt_param.TICKint)
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
                    ###############################################################
                    # Adjust Model Variables #
                    MDLsw_dwn = MDLsw_dwn[0:FXlen]; MDLnet_rad = MDLnet_rad[0:FXlen];
                    MDLtime = MDLtime[0:FXlen];
                    ###############################################################
                    SRADmin,SRADmax,SRADtick = atmo_ylim((MDLsw_dwn),plt_param.SRADrnd_bse,plt_param.TICKint)
                    NRADmin,NRADmax,NRADtick = atmo_ylim((MDLnet_rad),plt_param.SRADrnd_bse,plt_param.TICKint)
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
                    ###############################################################
                    # no model adjustment needed -- either equal to + (FXlen - 1) Fxcst or
                    # or less due to a download issue (i.e. MDLws = MDLws[0:FXlen])
                    ###############################################################
                    SRADmin,SRADmax,SRADtick = atmo_ylim((MDLsw_dwn),plt_param.SRADrnd_bse,plt_param.TICKint)
                    NRADmin,NRADmax,NRADtick = atmo_ylim((MDLnet_rad),plt_param.SRADrnd_bse,plt_param.TICKint)
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
                Shortwave Radiation Measurements
                ############################################################### """
                rad_fig, (sw_ax,net_ax,bias_ax) = plt.subplots(3,1,sharex = True, \
                          figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),facecolor = 'w',edgecolor ='w',clear='True')
                ###################################################################
                # Plot the Shortwave Downwelling Time Histories
                ###############################################
                sw_ax.plot(MDLtime,MDLtime*np.nan,'s-',color='gray',markersize=plt_param.PLTmsize,label='Obs (1Hr Avg)')
                sw_ax.plot(MDLtime,MDLtime*np.nan,'-',color='black',linewidth=plt_param.PLTlwidth,markersize=plt_param.PLTmsize,label='Obs')
                sw_ax.plot(MDLtime,MDLsw_dwn,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname+MDLstr)
                sw_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ####################################################################
                handles, labels = sw_ax.get_legend_handles_labels()
                order = [1,0,2]
                sw_ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], \
                               #loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_tri_ur,\
                               loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = (1.025,1.25),\
                               ncol = 3,frameon = False,columnspacing = plt_param.LEGcol_space);    
                ###################################################################
                sw_ax.set_ylim(SRADmin,SRADmax); sw_ax.set_yticks(SRADtick);
                sw_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                sw_ax.set_ylabel('Solar Radiation,\n W $\mathrm{m}^{-2}$',{'fontsize': plt_param.LABELfs});
                sw_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                ###################################################################
                #insert station info as title
                T = sw_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                T.set_position(plt_param.Tloc);
                             
                ###################################################################
                # Plot the Net Radiation Time Histories
                ###############################################
                net_ax.plot(MDLtime,MDLtime*np.nan,'s-',color='gray',markersize=plt_param.PLTmsize,label='Obs (1Hr Avg)')
                net_ax.plot(MDLtime,MDLtime*np.nan,'-',color='black',linewidth=plt_param.PLTlwidth,markersize=plt_param.PLTmsize,label='Obs')
                net_ax.plot(MDLtime,MDLnet_rad,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname+MDLstr)
                net_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ####################################################################
                net_ax.set_ylim(NRADmin,NRADmax); net_ax.set_yticks(NRADtick);
                net_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ###################################################################
                net_ax.set_ylabel('Net Radiation,\n W $\mathrm{m}^{-2}$',{'fontsize': plt_param.LABELfs});
                net_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                ###################################################################
                # 0/0 line
                net_ax.plot([Xmin,Xmax],[0,0],color='black',linewidth=.5)

                """ ###############################################################
                # Plot the Model Radiation Biases
                ############################################################### """
                # initate 0/0 line
                ###################################################################
                bias_ax.plot([Xmin,Xmax],[0,0],color='black',linewidth=1.5)
                bias_ax.plot(MDLtime,MDLtime*np.nan,'^-',color='forestgreen',markersize=plt_param.PLTmsize,\
                     linewidth=plt_param.PLTlwidth,label='Solar Radiation (RMSE: N/A)')
                bias_ax.plot(MDLtime,MDLtime*np.nan,'^-',color='dodgerblue',markersize=plt_param.PLTmsize,\
                     linewidth=plt_param.PLTlwidth,label='Net Radiation (RMSE: N/A)')
                bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                ###################################################################
                #bias_ax.legend(loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                #             ncol = 1,frameon = False);  
                #######################################################################
                bias_ax.set_xticks(Xtick); bias_ax.set_xticklabels(Xtick_lab,fontsize=plt_param.TICKfs);
                bias_ax.set_xlim(Xmin,Xmax)
                #######################################################################
                bias_ax.set_xlabel('Hour, UTC',{'fontsize': plt_param.LABELfs});
                bias_ax.xaxis.set_label_coords(plt_param.Xlab_pad_tri[0],plt_param.Xlab_pad_tri[1]);
                #######################################################################
                # Denote that Model Data was Not Availabel
                bias_ax.text(np.nanmean([Xmin,Xmax]),0,'Measurement Data Not Available',\
                             horizontalalignment='center',verticalalignment='bottom',color='red',\
                             fontweight='semibold',fontsize=plt_param.TITLEfs);
                #######################################################################
                #Define Bounds of Wind Speed - Model Bias to Determine Y-Bounds
                bias_ax.set_ylim(-100,100); bias_ax.set_yticks(np.arange(-100,100.1,50));
                bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                ####################################################################
                bias_ax.set_ylabel('Model Error,\n W $\mathrm{m}^{-2}$',{'fontsize': plt_param.LABELfs});    
                bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                ###################################################################
                bias_ax.legend(loc = plt_param.BIASpad,fontsize=plt_param.LEGfs,ncol = 2,\
                                      frameon = False);
                ###################################################################
                bias_ax.annotate('Mdl - Obs',xy=(0.13,0.345),xycoords='figure fraction',color='black',\
                             fontsize=plt_param.LEGfs,weight='bold');
                ###################################################################
                rad_fig.subplots_adjust(left = 0.125,right = 0.935,bottom = 0.155, top = 0.925)
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
                #######################################################################
                # add a logo
                imax = rad_fig.add_axes(plt_param.metLOGOaxcoords)
                imax.set_axis_off()
                imax.imshow(logo, aspect="equal")
                #######################################################################
                #Save the Developed Figure
                FILout = os.path.join(IMGpth,FILpfx+MDLname+'_SRAD_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                rad_fig.savefig(FILout,dpi=300); plt.close('all');
    
                """ Plot Extended Forecast Period if Applicable """
                ###################################################################
                if 'MDLtime_ext' in locals():
                    """ ###############################################################
                    Shortwave Radiation Measurements
                    ############################################################### """
                    rad_fig, (sw_ax,net_ax,bias_ax) = plt.subplots(3,1,sharex = True, \
                              figsize=(plt_param.FIGwdth,plt_param.FIGwdth*plt_param.FIGasp),facecolor = 'w',edgecolor ='w',clear='True')
                    ###################################################################
                    # Plot the Shortwave Downwelling Time Histories
                    ###############################################
                    sw_ax.plot(MDLtime_ext,MDLtime_ext*np.nan,'s-',color='gray',markersize=plt_param.PLTmsize,label='Obs (1Hr Avg)')
                    sw_ax.plot(MDLtime_ext,MDLtime_ext*np.nan,'-',color='black',linewidth=plt_param.PLTlwidth,markersize=plt_param.PLTmsize,label='Obs')
                    sw_ax.plot(MDLtime_ext,MDLsw_dwn_ext,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname+MDLstr)
                    sw_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                    ####################################################################
                    handles, labels = sw_ax.get_legend_handles_labels()
                    order = [1,0,2]
                    sw_ax.legend([handles[idx] for idx in order],[labels[idx] for idx in order], \
                                   #loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_tri_ur,\
                                   loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = (1.025,1.25),\
                                   ncol = 3,frameon = False,columnspacing = plt_param.LEGcol_space);    
                    ###################################################################
                    sw_ax.set_ylim(SRADmin_ext,SRADmax_ext); sw_ax.set_yticks(SRADtick_ext);
                    sw_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                    ###################################################################
                    sw_ax.set_ylabel('Solar Radiation,\n W $\mathrm{m}^{-2}$',{'fontsize': plt_param.LABELfs});
                    sw_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                    ###################################################################
                    #insert station info as title
                    T = sw_ax.set_title(MEAStitle,{'fontsize': plt_param.TITLEfs,'horizontalalignment':'left'},pad=2.5);
                    T.set_position(plt_param.Tloc);
                                 
                    ###################################################################
                    # Plot the Net Radiation Time Histories
                    ###############################################
                    net_ax.plot(MDLtime_ext,MDLtime_ext*np.nan,'s-',color='gray',markersize=plt_param.PLTmsize,label='Obs (1Hr Avg)')
                    net_ax.plot(MDLtime_ext,MDLtime_ext*np.nan,'-',color='black',linewidth=plt_param.PLTlwidth,markersize=plt_param.PLTmsize,label='Obs')
                    net_ax.plot(MDLtime_ext,MDLnet_rad_ext,'o-',color='red',markersize=plt_param.PLTmsize,linewidth=plt_param.PLTlwidth,label=MDLname+MDLstr)
                    net_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                    ####################################################################
                    net_ax.set_ylim(NRADmin,NRADmax); net_ax.set_yticks(NRADtick);
                    net_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                    ###################################################################
                    net_ax.set_ylabel('Net Radiation,\n W $\mathrm{m}^{-2}$',{'fontsize': plt_param.LABELfs});
                    net_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                    ###################################################################
                    # 0/0 line
                    net_ax.plot([Xmin_ext,Xmax_ext],[0,0],color='black',linewidth=.5)
                    ###################################################################
    
                    """ ###############################################################
                    # Plot the Model Radiation Biases
                    ############################################################### """
                    # initate 0/0 line
                    bias_ax.plot([Xmin_ext,Xmax_ext],[0,0],color='black',linewidth=1.5)
                    ###################################################################
                    bias_ax.plot(MDLtime_ext,MDLtime_ext*np.nan,'^-',color='forestgreen',markersize=plt_param.PLTmsize,\
                         linewidth=plt_param.PLTlwidth,label='Solar Radiation (RMSE: N/A)')
                    bias_ax.plot(MDLtime_ext,MDLtime_ext*np.nan,'^-',color='dodgerblue',markersize=plt_param.PLTmsize,\
                         linewidth=plt_param.PLTlwidth,label='Net Radiation (RMSE: N/A)')
                    bias_ax.grid('both',color='k',linestyle='-',linewidth=0.1,alpha=0.5)
                    ###################################################################
                    #bias_ax.legend(loc=1,fontsize=plt_param.LEGfs,bbox_to_anchor = plt_param.LEGpad_dual,\
                    #             ncol = 1,frameon = False);  
                    #######################################################################
                    bias_ax.set_xticks(Xtick_ext); bias_ax.set_xticklabels(Xtick_ext_lab,fontsize=plt_param.TICKfs);
                    bias_ax.set_xlim(Xmin_ext,Xmax_ext)
                    #######################################################################
                    bias_ax.set_xlabel('Hour, UTC',{'fontsize': plt_param.LABELfs});
                    bias_ax.xaxis.set_label_coords(plt_param.Xlab_pad_tri[0],plt_param.Xlab_pad_tri[1]);
                    #######################################################################
                    # Denote that Model Data was Not Availabel
                    bias_ax.text(np.nanmean([Xmin_ext,Xmax_ext]),0,'Measurement Data Not Available',\
                                 horizontalalignment='center',verticalalignment='bottom',color='red',\
                                 fontweight='semibold',fontsize=plt_param.TITLEfs);
                    #######################################################################
                    #Define Bounds of Wind Speed - Model Bias to Determine Y-Bounds
                    bias_ax.set_ylim(-100,100); bias_ax.set_yticks(np.arange(-100,100.1,50));
                    bias_ax.tick_params(axis="y", labelsize=plt_param.TICKfs)
                    ####################################################################
                    bias_ax.set_ylabel('Model Error,\n W $\mathrm{m}^{-2}$',{'fontsize': plt_param.LABELfs});    
                    bias_ax.yaxis.set_label_coords(plt_param.Ylab_pad[0],plt_param.Ylab_pad[1]);
                    ###################################################################
                    #bias_ax.legend(loc = plt_param.BIASpad,fontsize=plt_param.LEGfs,ncol = 2,\
                    #                      frameon = False);
                    ###################################################################
                    bias_ax.annotate('Mdl - Obs',xy=(0.13,0.345),xycoords='figure fraction',color='black',\
                                 fontsize=plt_param.LEGfs,weight='bold');
                    ###################################################################
                    rad_fig.subplots_adjust(left = 0.125,right = 0.935,bottom = 0.155, top = 0.925)
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
                    imax = rad_fig.add_axes(plt_param.metLOGOaxcoords)
                    imax.set_axis_off()
                    imax.imshow(logo, aspect="equal")
                    #######################################################################
                    #Save the Developed Figure
                    FILout = os.path.join(IMGpth_ext,FILpfx+MDLname+'_SRAD_'+Sid+'.'+DATE+'_INI'+INIpos[ini]+'.png')
                    rad_fig.savefig(FILout,dpi=300); plt.close('all');
                    ###############################################################
                    # Delete Extended Information to Preclude Incorrect Plotting
                    del MDLtime_ext
